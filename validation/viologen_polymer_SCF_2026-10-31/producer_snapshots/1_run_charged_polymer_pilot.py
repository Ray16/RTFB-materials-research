#!/usr/bin/env python3
"""Measure a route-bound ionic polymer SCF and initial force; qualify no G."""
import argparse
from contextlib import nullcontext
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
    def verify_sources():
        if any(hashlib.sha256(Path(p).read_bytes()).hexdigest()!=h for p,h in plan['source_sha256'].items()):
            raise ValueError('Charged-pilot source changed')
    verify_sources()
    if (plan['role']!='charged_polymer_SCF_gradient_resource_pilot' or plan['native_G_qualified']
            or plan['explicit_counterions_in_coordinates'] or len(os.sched_getaffinity(0))!=2
            or any(os.environ.get(k)!='2' for k in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','MKL_NUM_THREADS',
                'OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS','XTB_THREADS'))):raise ValueError('Declared ionic resource pilot and two-CPU execution required')
    record=dict(status='charged_polymer_SCF_gradient_pilot_running',started=time.time(),
        owner=dict(host=os.uname().nodename,pid=os.getpid(),process_start=Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()[19]),
        native_G_qualified=False,material_admission=False)
    atomic_write_json(case/'status.json',record)
    try:
        import numpy as np
        from ase.io import read
        from pyscf import gto
        from redox.qm import dft
        helper=load(sources/'helper.py','charged_pilot_corrected_screening')
        checkpoints=load(sources/'checkpoints.py','charged_pilot_SCF_checkpoints')
        model=load(ROOT/'scripts/pipeline/smooth_cds_control.py','charged_pilot_declared_CDS_model')
        monitor=load(ROOT/'scripts/pipeline/run_pi_stacking.py','charged_pilot_memory_monitor')
        if model.MODEL!=plan['solvent_model']:raise ValueError('The declared MeCN surface model differs')
        atoms=read(sources/'input.xyz')
        if len(atoms)!=plan['atoms'] or hashlib.sha256((sources/'input.xyz').read_bytes()).hexdigest()!=plan['selected_initial_geometry_sha256']:
            raise ValueError('Exact audited ionic geometry required')
        builder=checkpoints.recording_builder(dft._build_smd_mf,case/'SCF_evaluations',record['owner'])
        workspace=nullcontext()
        if plan['execution_profile']=='compact_28GiB':
            adapter=load(ROOT/'scripts/pipeline/df_occupied_workspace.py','charged_pilot_bounded_operators')
            workspace=adapter.compact_smd_execution(case/'execution',28)
        with helper.corrected_screening(sources/'kernel.cu',case/'screening_implementation'),monitor.MemoryMonitor() as memory,workspace,\
                model.smooth_cds(case/'model',ROOT/'results/validation/smooth_CDS_2026-10-25/pilot'):
            mol=gto.M(atom=list(zip(atoms.get_chemical_symbols(),atoms.positions.tolist())),unit='Angstrom',
                basis='def2-svp',charge=plan['charge'],spin=plan['mult']-1,verbose=0,max_memory=6000)
            mf=builder(mol,'r2scan','acetonitrile','gpu',disp='d4')
            mf.conv_tol=1e-12;mf.conv_tol_grad=1e-8;mf.max_cycle=200;mf.grids.level=3
            mf.callback=lambda e:atomic_write_json(case/'SCF_progress.json',dict(cycle=int(e['cycle'])+1,E_Ha=float(e['e_tot']),updated=time.time(),native_G_qualified=False))
            mf.kernel()
            if not mf.converged:raise ValueError('Initial charged-polymer SCF did not converge')
            spin=dft._orbital_diag(mf,plan['mult']-1)
            if abs(spin.get('spin_contam',math.inf))>.1:raise ValueError('Original spin gate failed')
            obj=mf.nuc_grad_method();obj.grid_response=True;value=obj.kernel()
            cpu=lambda x:np.asarray(x.get() if hasattr(x,'get') else x)
            gradient=cpu(value);parts={name:cpu(getattr(obj,'de_'+name)) for name in ('solute','solvent','cds')}
            if (gradient.shape!=(len(atoms),3) or not np.isfinite(gradient).all()
                    or np.max(abs(gradient-sum(parts.values())))>1e-10):raise ValueError('Finite full force components must close')
            np.save(case/'gradient.npy',gradient,allow_pickle=False)
            atomic_write_json(case/'evaluation.json',dict(coordinates_Bohr=mol.atom_coords().tolist(),
                charge=mol.charge,mult=mol.spin+1,nao=mol.nao,E_Ha=float(mf.e_tot),SCF_converged=True,spin=spin,
                gradient_Ha_Bohr=gradient.tolist(),gradient_components_Ha_Bohr={k:v.tolist() for k,v in parts.items()},
                gradient_max_Ha_Bohr=float(abs(gradient).max()),geometry_optimized=False,conformer_sampling_qualified=False,
                interaction_strength_computed=False,native_G_qualified=False))
        if not 0<memory.peak<=plan['maximum_peak_GPU_MiB']:raise ValueError('Charged pilot exceeds its prospective GPU envelope')
        verify_sources()
        record.update(status='charged_polymer_SCF_gradient_resource_pilot_finished_review_pending',
            observed_peak_GPU_MiB=memory.peak,gradient_max_Ha_Bohr=float(abs(gradient).max()),requires_solvent_optimization=True)
    except BaseException as error:
        record.update(status='charged_polymer_SCF_gradient_resource_pilot_failed_preserved',error=repr(error));raise
    finally:record['finished']=time.time();atomic_write_json(case/'status.json',record)


if __name__=='__main__':main()
