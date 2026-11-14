#!/usr/bin/env python
"""Two-model full quantum controls for an experimental smooth CDS surface.

Controls retain original force/Hessian thresholds. Results are method controls
in validation workspaces, never addressed primary candidate records.
"""
import argparse
import gc
from contextlib import nullcontext,redirect_stdout,redirect_stderr
import json
import os
from pathlib import Path
import time

import numpy as np

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mult',type=int,choices=(1,2),required=True)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--workspace-GiB',type=float,required=True)
    args=parser.parse_args()
    if args.workspace_GiB!=28:raise ValueError('Use the qualified 28 GiB workspace setting')
    from redox.workflows.resources import verify_gate_owner
    verify_gate_owner()
    caps=('OMP_NUM_THREADS','OMP_THREAD_LIMIT','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS','XTB_THREADS')
    if len(os.sched_getaffinity(0))!=2 or any(os.environ.get(k)!='2' for k in caps):
        raise ValueError('Two disjoint CPUs and native caps of two required')
    out=args.out.resolve()
    if not out.is_relative_to(ROOT/'results/validation'):
        raise ValueError('Experimental model controls must stay outside primary candidate directories')
    out.mkdir(parents=True,exist_ok=False)
    raw=ROOT/'calcs/validation'/out.parents[1].name/('spin'+str(args.mult));raw.mkdir(parents=True,exist_ok=False)
    adapter_path=ROOT/'scripts/pipeline/residual_response_solver.py'
    model_path=ROOT/'scripts/pipeline/smooth_cds_control.py'
    component_root=ROOT/'results/validation/smooth_CDS_2026-10-25/pilot'
    model=api.load_module('distinct_smooth_CDS_model',model_path)
    adapter=api.load_module('smooth_quantum_bounded_response',adapter_path)
    observer=api.load_module('smooth_control_memory',ROOT/'scripts/pipeline/run_pi_stacking.py')
    inputs=model.component_evidence(component_root)
    sources={str(p):api.digest(p) for p in (Path(__file__).resolve(),model_path,adapter_path,
        ROOT/'src/redox/qm/dft.py',ROOT/'src/redox/qm/cds_finite.py',ROOT/'scripts/validation/audit/native_partition.py')}
    xyz=raw/'input.xyz'
    xyz.write_text('4\nSinglet SMD model control\nC 0 0 0\nO 1.21 0 0\nH -.61 .94 0\nH -.61 -.94 0\n' if args.mult==1 else
        '4\nDoublet SMD model control\nC 0 0 0\nH 1.08 0 0\nH -.54 .935307436 0\nH -.54 -.935307436 0\n')
    report=dict(status='experimental_smooth_quantum_control_running',started=time.time(),mult=args.mult,
        owner=dict(host=os.uname().nodename,pid=os.getpid(),process_start=Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()[19]),
        source_sha256=sources,input_sha256=inputs,arms={},scientific_model_changed=True,
        primary_candidate_admission=False,raw_root=str(raw))
    atomic_write_json(out/'status.json',report)
    from ase.io import read
    from pyscf import gto
    from redox.qm import dft,qrrho
    from redox.qm.cds_finite import finer_cds_hessian
    from redox.workflows.geometry_checkpoints import checkpoint_callback
    equations=api.load_module('smooth_independent_nuclear_equations',ROOT/'scripts/validation/audit/native_partition.py')
    try:
        with observer.MemoryMonitor() as memory,(raw/'calculation.log').open('x') as log,\
                redirect_stdout(log),redirect_stderr(log):
            for arm in ('legacy','smooth'):
                folder=raw/arm;folder.mkdir()
                context=model.smooth_cds(folder/'model',component_root) if arm=='smooth' else nullcontext()
                with adapter.compact_smd_execution(folder/'execution',28),context:
                    with api.optimizer_settings(dict(max_opt_steps=100,geometry_SCF_conv_tol=1e-12),folder/'geometry_checkpoints'):
                        result=dft.dft_smd(xyz,0,args.mult,backend='gpu',do_opt=True,do_freq=False,
                            opt_out=folder/'opt.xyz',max_opt_steps=100,
                            conv_params=dict(convergence_gmax=1e-5,convergence_grms=6e-6,
                                convergence_energy=1e-8,convergence_dmax=.0005,convergence_drms=.0003,
                                callback=checkpoint_callback(folder/'geometry_checkpoints')))
                    if not all(result[k] for k in ('geometry_optimization_converged','optimization_SCF_converged','converged_smd','converged_gas')):
                        raise ValueError('Control geometry or electronic parent did not converge')
                    result.update(experimental_method_control=True,primary_candidate_admission=False,
                        CDS_profile=model.MODEL if arm=='smooth' else 'legacy_SMD_CDS')
                    atomic_write_json(folder/'quantum_result.json',result)
                    atoms=read(folder/'opt.xyz');symbols=atoms.get_chemical_symbols();coords=atoms.positions/.529177210903
                    counter=0
                    def solve(positions,gradient=True):
                        nonlocal counter
                        gc.collect();dft._free_gpu_pool()
                        dest=folder/'SCF_evaluations'/f'{counter:04d}';counter+=1;dest.mkdir(parents=True)
                        mol=gto.M(atom=list(zip(symbols,positions.tolist())),unit='Bohr',basis='def2-svp',
                            charge=0,spin=args.mult-1,verbose=0,max_memory=6000)
                        mf=dft._build_smd_mf(mol,'r2scan','acetonitrile','gpu',disp='d4')
                        mf.conv_tol=1e-12;mf.conv_tol_grad=1e-8;mf.max_cycle=200;mf.grids.level=3
                        mf.chkfile=str(dest/'SCF.h5');mf.kernel()
                        if not mf.converged:raise ValueError('Control SCF failed')
                        spin=dft._orbital_diag(mf,args.mult-1)
                        if abs(spin.get('spin_contam',float('inf')))>.1:
                            raise ValueError('Control spin contamination failed')
                        g=None
                        if gradient:
                            obj=mf.nuc_grad_method();obj.grid_response=True;v=obj.kernel()
                            g=np.asarray(v.get() if hasattr(v,'get') else v)
                        atomic_write_json(dest/'evaluation.json',dict(coordinates_Bohr=positions.tolist(),
                            E_Ha=float(mf.e_tot),gradient_Ha_Bohr=None if g is None else g.tolist(),
                            SCF_converged=True,spin=spin,arm=arm,qualified=False))
                        return mf,g
                    mf,gradient=solve(coords)
                    if abs(gradient).max()>.00045:raise ValueError('Smooth control is not stationary')
                    mf.cphf_grids=mf.grids.copy();hobj=mf.Hessian();hobj.grid_response=True;hobj.max_cycle=100
                    with finer_cds_hessian(1e-5):
                        value=hobj.kernel()
                    hessian=np.asarray(value.get() if hasattr(value,'get') else value).transpose(0,2,1,3).reshape(12,12)
                    np.save(folder/'hessian.npy',hessian);np.save(folder/'gradient.npy',gradient)
                    # All Cartesian full-energy forces, both original steps.
                    force_checks=[];column_checks=[]
                    selectors={(0,1),(1,2)}
                    for atom in range(4):
                        for axis in range(3):
                            for step in (.0005,.001):
                                plus,minus=coords.copy(),coords.copy();plus[atom,axis]+=step;minus[atom,axis]-=step
                                a,ga=solve(plus,(atom,axis) in selectors)
                                b,gb=solve(minus,(atom,axis) in selectors)
                                error=abs((a.e_tot-b.e_tot)/(2*step)-gradient[atom,axis])
                                force_checks.append(float(error))
                                if (atom,axis) in selectors:
                                    numerical=(ga-gb).ravel()/(2*step);analytic=hessian[:,3*atom+axis]
                                    column_checks.append(float(np.linalg.norm(numerical-analytic)/max(np.linalg.norm(analytic),1e-14)))
                                del a,b;dft._free_gpu_pool()
                    if max(force_checks)>1e-6 or max(column_checks)>.01:
                        raise ValueError('Full-energy force or Hessian-column validation failed')
                    nuclear=qrrho.nuclear_thermals(hessian,atoms.positions,atoms.numbers,1)
                    independent={str(k):equations.nuclear_equations(np.array(nuclear['frequencies_cm']),
                        atoms.positions,np.array(nuclear['atomic_masses_amu']),1,k) for k in (25,50,100)}
                    error=max(abs(independent[k]-nuclear['nuclear_by_cutoff_eV'][k]) for k in independent)
                    if error>1e-8:raise ValueError('Independent control nuclear partition failed')
                    atomic_write_json(folder/'nuclear.json',dict(nuclear,imposed_comparison_symmetry_number=1,
                        independent_nuclear_by_cutoff_eV=independent,maximum_reconstruction_error_eV=error,
                        primary_native_thermal_admission=False))
                    report['arms'][arm]=dict(E_smd_eV=result['e_smd_eV'],E_gas_eV=result['e_gas_eV'],
                        stationary_gradient_max_Ha_Bohr=float(abs(gradient).max()),
                        maximum_force_error_Ha_Bohr=max(force_checks),maximum_Hessian_column_relative_error=max(column_checks),
                        minimum_internal_frequency_cm=min(nuclear['frequencies_cm']),
                        nuclear_by_cutoff_eV=nuclear['nuclear_by_cutoff_eV'])
                    atomic_write_json(out/'status.json',report)
                    del mf,hobj;dft._free_gpu_pool()
        if not 0<memory.peak<=14336:raise ValueError('Small full control exceeded 14 GiB device bound')
        api.verify_sources(sources);api.verify_sources(inputs)
        report.update(status='two_model_quantum_force_Hessian_and_nuclear_controls_passed',
            observed_peak_GPU_MiB=memory.peak,output_sha256={str(p):api.digest(p) for p in raw.rglob('*')
                if p.is_file() and p.suffix in ('.json','.npy','.xyz','.py')})
    except BaseException as error:
        report.update(status='smooth_quantum_control_failed_preserved',error=repr(error));raise
    finally:
        report['finished']=time.time();atomic_write_json(out/'status.json',report)


if __name__=='__main__':
    main()
