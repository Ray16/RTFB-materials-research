#!/usr/bin/env python3
"""Optimize one exact ionic graph in MeCN with durable geometry evaluations."""
import argparse
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import time


def load(path,name):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--plan',type=Path,required=True);args=parser.parse_args()
    from redox.core.common import ROOT
    from redox.core.io import atomic_write_json
    from redox.workflows.resources import verify_gate_owner
    verify_gate_owner();case=args.plan.parent;sources=case.parent/'sources';plan=json.loads(args.plan.read_text())
    if plan['role']!='charged_polymer_solution_geometry_pilot' or plan['native_G_qualified']:
        raise ValueError('A separate route-bound charged geometry pilot is required')
    def check_sources():
        if any(hashlib.sha256(Path(p).read_bytes()).hexdigest()!=h for p,h in plan['source_sha256'].items()):
            raise ValueError('Geometry-pilot source changed')
    check_sources()
    if len(os.sched_getaffinity(0))!=2 or os.environ.get('NUMPY_MADVISE_HUGEPAGE')!='0' or any(
            os.environ.get(k)!='2' for k in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS','XTB_THREADS')):
        raise ValueError('Two disjoint CPUs and matching native caps required')
    record=dict(status='charged_polymer_solution_optimization_running',started=time.time(),
        owner=dict(host=os.uname().nodename,pid=os.getpid(),process_start=Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()[19]),
        native_G_qualified=False,material_admission=False)
    atomic_write_json(case/'status.json',record)
    try:
        import numpy as np
        from ase import Atoms
        from ase.io import read,write
        from pyscf import gto
        from pyscf.geomopt.geometric_solver import kernel
        from redox.qm import dft
        helper=load(sources/'helper.py','charged_geometry_screening')
        model=load(ROOT/'scripts/pipeline/smooth_cds_control.py','charged_geometry_CDS_model')
        adapter=load(ROOT/'scripts/pipeline/df_occupied_workspace.py','charged_geometry_memory_profile')
        monitor=load(ROOT/'scripts/pipeline/run_pi_stacking.py','charged_geometry_memory_monitor')
        if model.MODEL!=plan['solvent_model']:raise ValueError('Declared charged-polymer solvent model differs')
        atoms=read(sources/'input.xyz');symbols=atoms.get_chemical_symbols()
        cpu=lambda x:np.asarray(x.get() if hasattr(x,'get') else x)
        count=0
        def evaluation(data):
            nonlocal count
            count+=1;mol=data['mol'];gradient=cpu(data['gradients']);energy=float(data['energy'])
            if not np.isfinite(gradient).all() or not math.isfinite(energy):raise ValueError('Nonfinite geometry evaluation')
            point=case/'geometry_evaluations'/f'{count:04d}.json';point.parent.mkdir(parents=True,exist_ok=True)
            row=dict(owner=record['owner'],coordinates_Bohr=mol.atom_coords().tolist(),elements=symbols,
                charge=plan['charge'],mult=plan['mult'],E_Ha=energy,gradient_Ha_Bohr=gradient.tolist(),
                SCF_converged=bool(data['g_scanner'].converged),geometry_converged=False,qualified=False,updated=time.time())
            atomic_write_json(point,row);atomic_write_json(case/'latest_geometry.json',row)
        def build(mol,checkpoint):
            mf=dft._build_smd_mf(mol,'r2scan','acetonitrile','gpu',disp='d4')
            mf.conv_tol=1e-12;mf.conv_tol_grad=1e-8;mf.max_cycle=200;mf.grids.level=3;mf.chkfile=str(checkpoint)
            mf.callback=lambda e:atomic_write_json(case/'SCF_progress.json',dict(cycle=int(e['cycle'])+1,E_Ha=float(e['e_tot']),updated=time.time()))
            return mf
        with helper.corrected_screening(sources/'kernel.cu',case/'screening_implementation'),monitor.MemoryMonitor() as memory,\
                adapter.compact_smd_execution(case/'execution',28),\
                model.smooth_cds(case/'model',ROOT/'results/validation/smooth_CDS_2026-10-25/pilot'):
            mol=gto.M(atom=list(zip(symbols,atoms.positions.tolist())),unit='Angstrom',basis='def2-svp',
                charge=plan['charge'],spin=plan['mult']-1,verbose=0,max_memory=6000)
            mf=build(mol,case/'current_SCF.h5');gradient_method=mf.nuc_grad_method();gradient_method.grid_response=True
            converged,returned=kernel(gradient_method.as_scanner(),callback=evaluation,maxsteps=200,assert_convergence=True,
                convergence_energy=1e-6,convergence_grms=3e-4,convergence_gmax=4.5e-4,
                convergence_drms=1.2e-3,convergence_dmax=1.8e-3)
            atomic_write_json(case/'optimizer_return.json',dict(converged=bool(converged),coordinates_A=returned.atom_coords(unit='Angstrom').tolist(),
                elements=symbols,charge=plan['charge'],mult=plan['mult'],geometry_evaluations=count,qualified=False))
            final_atoms=Atoms(symbols,positions=returned.atom_coords(unit='Angstrom'));write(case/'returned.xyz',final_atoms)
            # A fresh solve at exactly the optimizer-returned coordinates avoids
            # attaching an earlier scanner energy to a different final geometry.
            final=build(returned,case/'final_SCF.h5');final.kernel()
            if not final.converged:raise ValueError('Final returned-geometry SCF did not converge')
            obj=final.nuc_grad_method();obj.grid_response=True;force=cpu(obj.kernel());spin=dft._orbital_diag(final,plan['mult']-1)
            if not np.isfinite(force).all() or abs(spin.get('spin_contam',math.inf))>.1:raise ValueError('Final force/spin gate failed')
            atomic_write_json(case/'final_evaluation.json',dict(coordinates_Bohr=returned.atom_coords().tolist(),E_Ha=float(final.e_tot),
                charge=plan['charge'],mult=plan['mult'],SCF_converged=True,spin=spin,gradient_Ha_Bohr=force.tolist(),
                gradient_max_Ha_Bohr=float(abs(force).max()),optimizer_converged=bool(converged),native_G_qualified=False))
        check_sources()
        if not 0<memory.peak<=31744:raise ValueError('Geometry pilot exceeds its proven 32-GB GPU profile')
        record.update(status='charged_polymer_solution_geometry_finished_independent_review_pending',
            optimizer_converged=bool(converged),maximum_final_gradient_Ha_Bohr=float(abs(force).max()),
            geometry_evaluations=count,observed_peak_GPU_MiB=memory.peak,Hessian_minimum_check_complete=False)
    except BaseException as error:
        record.update(status='charged_polymer_solution_optimization_failed_preserved',error=repr(error));raise
    finally:record['finished']=time.time();atomic_write_json(case/'status.json',record)


if __name__=='__main__':main()
