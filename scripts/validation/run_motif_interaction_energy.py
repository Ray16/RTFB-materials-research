#!/usr/bin/env python3
"""Calculate frozen A/B/AB energies with unchanged r2SCAN-D4/SVP SCF equations."""
import argparse
import gc
import importlib.util
import json
import math
import os
from pathlib import Path
import time


def load(path,name):
    definition=importlib.util.spec_from_file_location(name,path);module=importlib.util.module_from_spec(definition);definition.loader.exec_module(module);return module


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);a=p.parse_args()
    from redox.core.common import ROOT
    from redox.core.io import atomic_write_json
    from redox.workflows.resources import verify_gate_owner
    verify_gate_owner();case=a.plan.parent;sources=case.parent/'sources';plan=json.loads(a.plan.read_text())
    sha=load(sources/'hashing.py','motif_frozen_hashes').sha256_file
    def verify():
        if any(sha(path)!=h for path,h in plan['source_sha256'].items()):raise ValueError('Frozen motif calculation inputs changed')
    if plan['role']!='frozen_charged_motif_energy' or plan['native_G_qualified'] or plan['geometry_optimized']:
        raise ValueError('Exact electronic-only frozen charged motif control required')
    verify();owner=dict(host=os.uname().nodename.split('.')[0],pid=os.getpid(),process_start=Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()[19])
    if len(os.sched_getaffinity(0))!=2 or any(os.environ.get(k)!='2' for k in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS','XTB_THREADS')):
        raise ValueError('Two disjoint CPUs and consistent native thread caps required')
    stage=dict(status='frozen_charged_motif_energy_running',owner=owner,started=time.time(),native_G_qualified=False)
    atomic_write_json(case/'status.json',stage)
    try:
        from ase.io import read
        from pyscf import gto
        from redox.qm import dft
        helper=load(sources/'helper.py','motif_corrected_screening');checkpoints=load(sources/'checkpoints.py','motif_full_SCF_checkpoints')
        model=load(ROOT/'scripts/pipeline/smooth_cds_control.py','motif_matching_smooth_model')
        monitor=load(ROOT/'scripts/pipeline/run_pi_stacking.py','motif_GPU_monitor')
        if model.MODEL!=plan['solvent_model']:raise ValueError('Parent-model context changed')
        energies={};method=plan['method']
        with helper.corrected_screening(sources/'kernel.cu',case/'screening_implementation'),monitor.MemoryMonitor() as memory,\
                model.smooth_cds(case/'model',ROOT/'results/validation/smooth_CDS_2026-10-25/pilot'):
            for task in plan['tasks']:
                label=task['id'];atoms=read(task['geometry']);energies[label]={}
                for medium in plan['media']:
                    mol=gto.M(atom=list(zip(atoms.get_chemical_symbols(),atoms.positions.tolist())),basis=method['basis'],
                        charge=task['charge'],spin=task['mult']-1,verbose=0,max_memory=6000)
                    def build():
                        if medium=='acetonitrile':return dft._build_smd_mf(mol,method['xc'],medium,'gpu',disp=method['disp'])
                        mf=dft._dft_module('gpu').RKS(mol).density_fit();mf.xc=method['xc'];mf.disp=method['disp'];return mf
                    mf=checkpoints.recording_builder(build,case/'SCF_evaluations'/label/medium,owner,save_density=False)()
                    mf.conv_tol=method['conv_tol'];mf.conv_tol_grad=method['conv_tol_grad'];mf.max_cycle=method['max_cycle'];mf.grids.level=method['grid_level']
                    mf.kernel()
                    if not mf.converged or not math.isfinite(float(mf.e_tot)):raise ValueError('Frozen motif SCF did not converge')
                    spin=dft._orbital_diag(mf,0)
                    if abs(spin.get('spin_contam',float('inf')))>.1:raise ValueError('Frozen motif spin gate failed')
                    energies[label][medium]=dict(E_Ha=float(mf.e_tot),SCF_converged=True,spin=spin,
                        charge=task['charge'],mult=task['mult'],geometry_sha256=sha(task['geometry']))
                    atomic_write_json(case/'energies_partial.json',dict(energies=energies,native_G_qualified=False))
                    del mf,mol;gc.collect();dft._free_gpu_pool()
        if not 0<memory.peak<=31744:raise ValueError('Frozen motif resource pilot exceeded its GPU bound')
        verify();atomic_write_json(case/'result.json',dict(energies=energies,method=method,solvent_model=plan['solvent_model'],
            interaction_energy_qualified=False,native_G_qualified=False,basis_superposition_corrected=False,electronic_coupling_computed=False))
        stage.update(status='frozen_charged_motif_energy_finished_independent_review_pending',observed_peak_GPU_MiB=memory.peak)
    except BaseException as error:stage.update(status='frozen_charged_motif_energy_failed_preserved',error=repr(error));raise
    finally:stage['finished']=time.time();atomic_write_json(case/'status.json',stage)


if __name__=='__main__':main()
