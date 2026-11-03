#!/usr/bin/env python3
"""One source-frozen SCF/gradient on a prospective small GPU; no Hessian or G."""
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import time


def load(path,name):
    s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--plan',type=Path,required=True);args=parser.parse_args()
    from redox.core.common import ROOT
    from redox.core.io import atomic_write_json
    from redox.workflows.resources import verify_gate_owner
    verify_gate_owner();case=args.plan.parent;sources=case.parent/'sources';plan=json.loads(args.plan.read_text())
    profile=load(sources/'resource_profile.py','fixed_gradient_only_resource_profile');profile.validate_plan(plan)
    caps=('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','XTB_THREADS')
    if len(os.sched_getaffinity(0))!=2 or any(os.environ.get(k)!='2' for k in caps):
        raise ValueError('Two allocated CPUs and all six native caps=2 required')
    def verify():
        if any(hashlib.sha256(Path(p).read_bytes()).hexdigest()!=h for p,h in plan['source_sha256'].items()):
            raise ValueError('Fixed-gradient source snapshot changed')
    verify();owner=dict(host=os.uname().nodename.split('.')[0],pid=os.getpid(),
        process_start=Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()[19])
    stage=dict(status='counterion_fixed_gradient_running',owner=owner,started=time.time(),native_G_qualified=False)
    atomic_write_json(case/'status.json',stage)
    try:
        import numpy as np
        from ase.io import read
        from pyscf import gto
        from redox.qm import dft
        atoms=read(sources/'input.xyz')
        if len(atoms)!=42:raise ValueError('Exact reference atom count required')
        checkpoints=load(sources/'checkpoints.py','fixed_gradient_actual_SCF_checkpoints')
        screening=load(sources/'helper.py','fixed_gradient_corrected_screening')
        monitor=load(ROOT/'scripts/pipeline/run_pi_stacking.py','fixed_gradient_memory_monitor')
        builder=checkpoints.recording_builder(dft._build_smd_mf,case/'SCF_evaluations',owner,save_density=False)
        with screening.corrected_screening(sources/'kernel.cu',case/'screening_implementation'),monitor.MemoryMonitor() as memory:
            mol=gto.M(atom=list(zip(atoms.get_chemical_symbols(),atoms.positions.tolist())),unit='Angstrom',
                basis='def2-svp',charge=1,spin=0,verbose=0,max_memory=6000)
            if mol.nao_nr()!=412:raise ValueError('Exact primary geometry-basis dimension required')
            mf=builder(mol,'r2scan','acetonitrile','gpu',disp='d4')
            mf.conv_tol=1e-12;mf.conv_tol_grad=1e-8;mf.max_cycle=200;mf.grids.level=3
            mf.callback=lambda e:atomic_write_json(case/'progress.json',dict(cycle=int(e['cycle'])+1,E_Ha=float(e['e_tot']),updated=time.time()))
            mf.kernel()
            if not mf.converged:raise ValueError('Fixed-gradient SCF did not converge')
            derivative=mf.nuc_grad_method();derivative.grid_response=True;value=derivative.kernel()
            gradient=np.asarray(value.get() if hasattr(value,'get') else value)
            spin=dft._orbital_diag(mf,0)
            if not np.isfinite(gradient).all() or abs(spin.get('spin_contam',float('inf')))>.1:
                raise ValueError('Finite gradient and unchanged spin gate required')
            atomic_write_json(case/'evaluation.json',dict(E_Ha=float(mf.e_tot),gradient_Ha_Bohr=gradient.tolist(),
                coordinates_Bohr=mol.atom_coords().tolist(),SCF_converged=True,spin=spin,
                Hessian_computed=False,native_G_qualified=False))
        if not 0<memory.peak<=plan['maximum_peak_GPU_MiB']:raise ValueError('Prospective T4 memory headroom failed')
        verify();stage.update(status='counterion_fixed_gradient_finished_review_pending',observed_peak_GPU_MiB=memory.peak)
    except BaseException as error:
        stage.update(status='counterion_fixed_gradient_failed_preserved',error=repr(error));raise
    finally:
        stage['finished']=time.time();atomic_write_json(case/'status.json',stage)


if __name__=='__main__':main()
