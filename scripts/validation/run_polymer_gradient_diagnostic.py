#!/usr/bin/env python3
"""Repeat three failed-column SCFs and retain full gradient components; no G."""
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import time

import numpy as np


def load(path,name):
    s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);a=p.parse_args()
    from redox.core.common import ROOT
    sha256_file=lambda path:hashlib.sha256(Path(path).read_bytes()).hexdigest()
    from redox.core.io import atomic_write_json
    from redox.workflows.resources import verify_gate_owner,host_capacity,gpu_suitability
    from pyscf import gto
    from redox.qm import dft
    verify_gate_owner();plan=json.loads(a.plan.read_text());case=a.plan.parent
    if plan['role']!='polymer_gradient_diagnostic' or any(sha256_file(k)!=v for k,v in plan['source_sha256'].items()):raise ValueError('Diagnostic source identity changed')
    if len(os.sched_getaffinity(0))!=2 or any(os.environ.get(k)!='2' for k in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS','XTB_THREADS')):raise ValueError('Two CPUs and matching thread caps required')
    host_capacity(2,minimum_ram_gib=112);gpu_suitability(int(os.environ['CUDA_VISIBLE_DEVICES']))
    owner=dict(host=os.uname().nodename,pid=os.getpid(),process_start=Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()[19])
    record=dict(status='polymer_gradient_diagnostic_running',owner=owner,started=time.time(),native_G_qualified=False,SCF_points=[])
    atomic_write_json(case/'status.json',record)
    # Gradient-only work does not exercise a Hessian response solve. Use the
    # identical underlying execution operators without its response-use gate.
    adapter=load(ROOT/'scripts/pipeline/df_occupied_workspace.py','diagnostic_original_workspace_adapter')
    model=load(ROOT/'scripts/pipeline/smooth_cds_control.py','diagnostic_original_smooth_model')
    monitor=load(ROOT/'scripts/pipeline/run_pi_stacking.py','diagnostic_original_memory_monitor')
    checkpoint=load(case.parent/'sources/checkpoints.py','diagnostic_durable_checkpoints')
    dft._build_smd_mf=checkpoint.recording_builder(dft._build_smd_mf,case/'SCF_evaluations',owner)
    def cpu(x):return np.asarray(x.get() if hasattr(x,'get') else x)
    try:
        with monitor.MemoryMonitor() as memory,adapter.compact_smd_execution(case/'execution',28),\
                model.smooth_cds(case/'model',ROOT/'results/validation/smooth_CDS_2026-10-25/pilot'):
            for index,point in enumerate(plan['points']):
                folder=case/f'point{index}';folder.mkdir();original=json.loads(Path(point['evaluation']).read_text())
                mol=gto.M(atom=list(zip(plan['symbols'],original['coordinates_Bohr'])),unit='Bohr',basis='def2-svp',charge=0,spin=0,verbose=0,max_memory=6000)
                mf=dft._build_smd_mf(mol,'r2scan','acetonitrile','gpu',disp='d4')
                mf.conv_tol=1e-12;mf.conv_tol_grad=1e-8;mf.max_cycle=200;mf.grids.level=3
                mf.callback=lambda e:atomic_write_json(folder/'progress.json',dict(cycle=int(e['cycle'])+1,E_Ha=float(e['e_tot']),updated=time.time(),native_G_qualified=False))
                mf.kernel()
                if not mf.converged:raise ValueError('Diagnostic SCF did not converge')
                obj=mf.nuc_grad_method();obj.grid_response=True;gradient=cpu(obj.kernel())
                parts={name:cpu(getattr(obj,'de_'+name)) for name in ('solute','solvent','cds')}
                if (not np.isfinite(gradient).all() or not all(np.isfinite(v).all() for v in parts.values())
                        or np.max(abs(gradient-sum(parts.values())))>1e-10):raise ValueError('Gradient components do not close')
                np.save(folder/'gradient.npy',gradient,allow_pickle=False)
                for name,value in parts.items():np.save(folder/(name+'.npy'),value,allow_pickle=False)
                row=dict(source_point=point['source_point'],E_Ha=float(mf.e_tot),source_E_Ha=original['E_Ha'],
                    energy_difference_Ha=float(mf.e_tot)-original['E_Ha'],
                    maximum_gradient_difference_Ha_Bohr=float(np.max(abs(gradient-np.asarray(original['gradient_Ha_Bohr'])))),
                    completed=time.time(),gradient_grid_response=True,SCF_converged=True)
                atomic_write_json(folder/'result.json',row);record['SCF_points'].append(row);atomic_write_json(case/'status.json',record)
                del obj,mf;import gc;gc.collect();dft._free_gpu_pool()
        if not 0<memory.peak<=31744:raise ValueError('Diagnostic GPU ceiling failed')
        record.update(status='polymer_gradient_components_written_independent_analysis_pending',observed_peak_GPU_MiB=memory.peak)
    except BaseException as error:record.update(status='polymer_gradient_diagnostic_failed_preserved',error=repr(error));raise
    finally:record['finished']=time.time();atomic_write_json(case/'status.json',record)


if __name__=='__main__':main()
