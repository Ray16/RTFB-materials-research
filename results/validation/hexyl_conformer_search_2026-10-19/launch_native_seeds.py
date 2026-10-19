"""Use all free suitable local GPUs up to the two predeclared seed calculations."""
import json
import os
from pathlib import Path
import socket
import subprocess
import sys
import time

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows.resources import allocation_guard,cpu_allocation,free_slots,gpu_suitability,host_capacity

base=Path(__file__).resolve().parent;host=socket.gethostname().split('.')[0]
plan=json.loads((base/'native_seed_plan.json').read_text())
preflight=json.loads((base/f'native_preflight_{host}.json').read_text())
assert preflight['ok'] and preflight['python']==sys.executable
receipt=base/'native_launch.json';assert not receipt.exists()
with allocation_guard():
    suitable=[]
    for _,gpu in free_slots([host]):
        try:device=gpu_suitability(gpu)
        except RuntimeError:continue
        suitable.append((gpu,device))
    suitable=suitable[:len(plan['tasks'])];assert suitable
    resources=host_capacity(2*len(suitable),minimum_ram_gib=32*len(suitable))
    cpus=cpu_allocation(2*len(suitable));report=dict(host=host,resources=resources,
        workers=[],plan=str(base/'native_seed_plan.json'),status='submitting',created=time.time())
    env=dict(os.environ,PYTHONPATH=str(ROOT/'src'),REDOX_HESSIAN_MEMORY_GIB='24')
    for key in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS','XTB_THREADS'):env[key]='2'
    for index,(gpu,device) in enumerate(suitable):
        assigned=cpus[2*index:2*index+2]
        reservation=base/f'native_reserved_{host}_gpu{gpu}'
        argv=[sys.executable,'-u','-m','redox.workflows.gpu','--gpu',str(gpu),
            '--cpus',','.join(map(str,assigned)),'--threads','2','--out',str(reservation),
            '--',sys.executable,'-u',str(base/'native_seed_worker.py'),'--task',str(index)]
        with (base/f'native_controller_{host}_gpu{gpu}.log').open('w') as stream:
            proc=subprocess.Popen(argv,cwd=ROOT,env=env,stdin=subprocess.DEVNULL,
                stdout=stream,stderr=subprocess.STDOUT,start_new_session=True,
                preexec_fn=lambda assigned=assigned:os.sched_setaffinity(0,assigned))
        report['workers'].append(dict(key=plan['tasks'][index]['key'],gpu=gpu,device=device,
            cpu_affinity=assigned,controller_pid=proc.pid,
            process_start=Path(f'/proc/{proc.pid}/stat').read_text().rsplit(')',1)[1].split()[19],
            argv=argv,reservation=str(reservation),submitted=time.time()))
        atomic_write_json(receipt,report)
        time.sleep(.5)
    report['status']='submitted_predictor_verification_pending';atomic_write_json(receipt,report)
print('Submitted',len(suitable),'private hexyl seed calculations')
