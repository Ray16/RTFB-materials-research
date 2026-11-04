#!/usr/bin/env python3
"""Compare NumPy huge-page advice on identical saved data without starting QM.

Run fresh, one-CPU processes in the sequence off/on/off. A timed-out sample is
censored performance evidence, not a numerical failure or a scientific result.
No shared kernel setting or running worker is modified.
"""
import argparse
import faulthandler
import hashlib
import json
import os
from pathlib import Path
import resource
import subprocess
import sys
import time

from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json


def worker(raw,output,flag):
    faulthandler.dump_traceback_later(5,repeat=True)
    import numpy as np
    import h5py
    start=time.monotonic();usage=resource.getrusage(resource.RUSAGE_SELF)
    identity=json.loads((raw/'execution/rks_make_h1_state_identity.json').read_text())
    atomic_write_json(output/'started.json',dict(hugepage_advice=flag,numpy=np.__version__,h5py=h5py.__version__,HDF5=h5py.version.hdf5_version,CPU_affinity=sorted(os.sched_getaffinity(0)),pid=os.getpid()))
    atomic_write_json(output/'progress.json',dict(phase='read_saved_orbitals',qualified_scientific_result=False))
    with h5py.File(raw/'SCF_evaluations/0000/SCF.h5','r') as handle:coeff=handle['scf/mo_coeff'][()]
    digest=hashlib.sha256(coeff.tobytes()).hexdigest()
    if coeff.dtype!=np.float64 or digest!=identity['mo_coeff']['sha256']:raise ValueError('Control does not reproduce the source orbital payload')
    atomic_write_json(output/'progress.json',dict(phase='scan_saved_first_derivatives',orbital_payload_sha256=digest,qualified_scientific_result=False))
    h1=np.load(raw/'execution/rks_make_h1_0.npy',mmap_mode='r',allow_pickle=False)
    finite=bool(np.isfinite(h1).all());final=resource.getrusage(resource.RUSAGE_SELF)
    if not finite:raise ValueError('Nonfinite source first-derivative tensor')
    result=dict(status='allocation_control_completed',hugepage_advice=flag,wall_seconds=time.monotonic()-start,
        user_seconds=final.ru_utime-usage.ru_utime,system_seconds=final.ru_stime-usage.ru_stime,
        peak_RSS_bytes=final.ru_maxrss*1024,orbital_payload_sha256=digest,H1_shape=list(h1.shape),all_finite=finite,native_G_qualified=False)
    atomic_write_json(output/'result.json',result);faulthandler.cancel_dump_traceback_later()


def compare(raw,output,seconds):
    from redox.workflows.resources import allocation_guard,host_capacity,cpu_allocation
    with allocation_guard():host_capacity(1,minimum_ram_gib=4);os.sched_setaffinity(0,cpu_allocation(1))
    if output.exists():raise FileExistsError('Preserve earlier allocation comparisons')
    output.mkdir(parents=True)
    paths=[raw/'SCF_evaluations/0000/SCF.h5',raw/'execution/rks_make_h1_0.npy',raw/'execution/rks_make_h1_state_identity.json',Path(__file__).resolve()]
    hashes={str(p):sha256_file(p) for p in paths};rows=[]
    for index,flag in enumerate([0,1,0]):
        folder=output/f'{index:02d}_hugepage_{flag}';folder.mkdir()
        env=dict(os.environ,NUMPY_MADVISE_HUGEPAGE=str(flag))
        for key in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','XTB_THREADS'):env[key]='1'
        command=[sys.executable,str(Path(__file__).resolve()),'--worker','--flag',str(flag),'--raw',str(raw),'--out',str(folder)]
        started=time.monotonic()
        with (folder/'process.log').open('x') as log:
            try:
                run=subprocess.run(command,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=seconds,check=False)
                row=dict(flag=flag,status='completed' if run.returncode==0 else 'process_failed',returncode=run.returncode)
            except subprocess.TimeoutExpired:
                row=dict(flag=flag,status='observation_budget_expired',observation_seconds=seconds)
        row['elapsed_seconds']=time.monotonic()-started
        result=folder/'result.json'
        if result.exists():row['result']=json.loads(result.read_text())
        rows.append(row);atomic_write_json(output/'progress.json',dict(samples=rows,complete=False))
    if any(sha256_file(p)!=h for p,h in hashes.items()):raise ValueError('Control source changed')
    completed=[row['result'] for row in rows if row['status']=='completed']
    if len({row['orbital_payload_sha256'] for row in completed})!=1:raise ValueError('Completed allocation controls read different payloads')
    result=dict(status='source_bound_numpy_allocation_comparison_completed',source_sha256=hashes,sequence=[0,1,0],samples=rows,
        disabled_samples_all_completed=all(row['status']=='completed' for row in rows if row['flag']==0),
        enabled_observation_censored=any(row['status']=='observation_budget_expired' for row in rows if row['flag']==1),
        numerical_QM_protocol_changed=False,shared_kernel_policy_changed=False,original_scientific_workers_signalled=False,native_G_qualified=False,
        references={'NumPy':'https://numpy.org/doc/2.0/reference/global_state.html','kernel':'https://www.kernel.org/doc/html/v6.7/admin-guide/mm/transhuge.html'})
    atomic_write_json(output/'comparison.json',result);return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--raw',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--seconds',type=float,default=20);p.add_argument('--worker',action='store_true');p.add_argument('--flag',type=int,choices=(0,1));a=p.parse_args()
    if not 5<=a.seconds<=60:p.error('Use a finite 5--60 s per-sample observation budget')
    if a.worker:
        if a.flag is None or os.environ.get('NUMPY_MADVISE_HUGEPAGE')!=str(a.flag):p.error('Explicit import-time allocation flag required')
        worker(a.raw.resolve(),a.out.resolve(),a.flag)
    else:
        r=compare(a.raw.resolve(),a.out.resolve(),a.seconds);print(json.dumps(dict(status=r['status'],samples=[dict(flag=x['flag'],status=x['status'],seconds=x.get('result',{}).get('wall_seconds')) for x in r['samples']])))
