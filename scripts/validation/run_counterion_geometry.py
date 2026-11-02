#!/usr/bin/env python3
"""Optimize one explicitly defined chloride-contact hypothesis in primary SMD.

No smooth-surface polymer model is substituted. Every geometry/SCF evaluation
and optimizer return is retained; only fresh final-SCF/gradient checks can
qualify a geometry. Native G, association and populations remain pending.
"""
import argparse
import gc
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import time


def load(path,name):
    definition=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(definition);definition.loader.exec_module(module);return module


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);a=p.parse_args()
    from redox.core.common import ROOT
    from redox.core.io import atomic_write_json
    from redox.workflows.resources import verify_gate_owner
    verify_gate_owner();case=a.plan.parent;sources=case.parent/'sources';plan=json.loads(a.plan.read_text())
    if (plan['role']!='primary_counterion_geometry_control' or (plan['charge'],plan['mult'])!=(1,1)
            or plan['parent_formal_core_charge']!=2 or plan['counterion_charge']!=-1
            or plan['native_G_qualified'] or plan['material_admission']):
        raise ValueError('Exact oxidized-core chloride sensitivity scope required')
    def verify():
        if any(hashlib.sha256(Path(p).read_bytes()).hexdigest()!=h for p,h in plan['source_sha256'].items()):
            raise ValueError('Counterion geometry source snapshot changed')
    verify();owner=dict(host=os.uname().nodename.split('.')[0],pid=os.getpid(),
        process_start=Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()[19])
    stage=dict(status='counterion_geometry_control_running',owner=owner,started=time.time(),native_G_qualified=False)
    atomic_write_json(case/'status.json',stage)
    try:
        import numpy as np
        from ase.io import read
        from pyscf import gto
        from redox.qm import dft
        from redox.validation.protonation import hydrogen_topology_status
        checkpoints=load(sources/'checkpoints.py','counterion_full_checkpoints')
        screening=load(sources/'helper.py','counterion_corrected_screening')
        monitor=load(ROOT/'scripts/pipeline/run_pi_stacking.py','counterion_GPU_monitor')
        builder=dft._build_smd_mf;optimizer=dft._geomopt_fn
        def strict(*args,**kwargs):
            mf=builder(*args,**kwargs);mf.conv_tol=1e-12;mf.conv_tol_grad=1e-8;mf.max_cycle=200;mf.grids.level=3
            return mf
        dft._build_smd_mf=checkpoints.recording_builder(strict,case/'SCF_evaluations',owner,save_density=False)
        dft._geomopt_fn=lambda backend:checkpoints.recording_optimizer(optimizer(backend),case/'optimizer_returns',owner)
        try:
            with screening.corrected_screening(sources/'kernel.cu',case/'screening_implementation'),monitor.MemoryMonitor() as memory:
                with (case/'optimization.log').open('x') as log:
                    from contextlib import redirect_stdout,redirect_stderr
                    with redirect_stdout(log),redirect_stderr(log):
                        result=dft.dft_smd(case/'input.xyz',1,1,backend='gpu',solvent='acetonitrile',
                            do_opt=True,do_freq=False,do_smd=False,do_gas=False,opt_out=case/'optimizer_return.xyz',
                            max_opt_steps=400,assert_convergence=False,
                            conv_params=plan.get('optimizer_convergence'))
                atomic_write_json(case/'optimizer_result.json',result)
                if not result['optimized']:raise ValueError('Optimizer did not converge; all returned coordinates remain unqualified')
                if hydrogen_topology_status(plan['smiles'],case/'optimizer_return.xyz')!='ok':
                    raise ValueError('Counterion optimization changed the exact intended ionic/covalent graph')
                gc.collect();dft._free_gpu_pool()
                atoms=read(case/'optimizer_return.xyz');mol=gto.M(atom=list(zip(atoms.get_chemical_symbols(),atoms.positions.tolist())),
                    basis='def2-svp',charge=1,spin=0,verbose=0,max_memory=6000)
                mf=dft._build_smd_mf(mol,'r2scan','acetonitrile','gpu',disp='d4');mf.kernel()
                if not mf.converged:raise ValueError('Fresh optimizer-return SCF did not converge')
                obj=mf.nuc_grad_method();obj.grid_response=True;value=obj.kernel()
                gradient=np.asarray(value.get() if hasattr(value,'get') else value);spin=dft._orbital_diag(mf,0)
                atomic_write_json(case/'final_evaluation.json',dict(coordinates_Bohr=mol.atom_coords().tolist(),
                    elements=atoms.get_chemical_symbols(),E_Ha=float(mf.e_tot),gradient_Ha_Bohr=gradient.tolist(),
                    spin=spin,charge=1,mult=1,SCF_converged=True,method=plan['method'],geometry_qualified=False,native_G_qualified=False))
                if not np.isfinite(gradient).all() or abs(gradient).max()>.00045 or abs(spin.get('spin_contam',float('inf')))>.1:
                    raise ValueError('Fresh final stationarity/spin gate failed')
            if not 0<memory.peak<=31744:raise ValueError('Counterion geometry GPU memory bound exceeded')
        finally:dft._build_smd_mf=builder;dft._geomopt_fn=optimizer
        verify();stage.update(status='counterion_geometry_control_finished_independent_review_pending',
            observed_peak_GPU_MiB=memory.peak,maximum_gradient_Ha_Bohr=float(abs(gradient).max()))
    except BaseException as error:
        stage.update(status='counterion_geometry_control_failed_preserved',error=repr(error));raise
    finally:
        stage['finished']=time.time();atomic_write_json(case/'status.json',stage)


if __name__=='__main__':main()
