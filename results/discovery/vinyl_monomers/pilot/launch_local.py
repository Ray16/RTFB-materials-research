"""Add all free suitable local GPUs to the shared, bounded molecular queue."""
import json
import os
from pathlib import Path
import socket
import site
import subprocess
import sys
import time

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows.resources import allocation_guard, cpu_allocation, free_slots, gpu_suitability, host_capacity

base=Path(__file__).resolve().parent
host=socket.gethostname().split('.')[0]
plan=json.loads((base/'plan.json').read_text())
preflight=json.loads((base/f'{host}_preflight.json').read_text())
assert preflight['ok'] and preflight['python']==sys.executable
receipt=base/f'launch_{host}.json'
assert not receipt.exists(), 'Existing launch is preserved; inspect ownership before any resume'
assert len(plan['tasks'])==9
with allocation_guard():
    slots=free_slots([host])
    suitable=[]
    for _,gpu in slots:
        try:device=gpu_suitability(gpu)
        except RuntimeError:continue
        if device['memory_MiB']>=80*1024:suitable.append((gpu,device))
    suitable=suitable[:len(plan['tasks'])]
    assert suitable, 'No currently free suitable 32-GB-class devices'
    resources=host_capacity(2*len(suitable),minimum_ram_gib=32*len(suitable))
    cpus=cpu_allocation(2*len(suitable))
    env=dict(os.environ,PYTHONPATH=str(ROOT/'src'),REDOX_HESSIAN_MEMORY_GIB='24',
             FAIRCHEM_CACHE_DIR='/nfs/lambda_stor_01/homes/rzhu/.cache/fairchem')
    env['PATH']=str(Path(sys.executable).parent)+os.pathsep+env.get('PATH','')
    env['LD_LIBRARY_PATH']=':'.join(str(p) for b in site.getsitepackages() for p in Path(b).glob('nvidia/*/lib'))+':'+env.get('LD_LIBRARY_PATH','')
    for key in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS','XTB_THREADS'):
        env[key]='2'
    report=dict(status='submitting_shared_queue_workers',host=host,resources=resources,
        workers=[],source_plan=str(base/'plan.json'),created=time.time())
    atomic_write_json(receipt,report)
    for i,(gpu,device) in enumerate(suitable):
        assigned=cpus[2*i:2*i+2]
        argv=[sys.executable,'-u','-m','redox.workflows.gpu','--gpu',str(gpu),
              '--cpus',','.join(map(str,assigned)),'--threads','2',
              '--out',str(base/f'reserved_{host}_gpu{gpu}'),'--',
              sys.executable,'-u',str(base/'queue_worker.py')]
        log=base/f'controller_{host}_gpu{gpu}.log'
        with log.open('w') as stream:
            proc=subprocess.Popen(argv,cwd=ROOT,env=env,stdin=subprocess.DEVNULL,
                stdout=stream,stderr=subprocess.STDOUT,start_new_session=True,
                preexec_fn=lambda assigned=assigned:os.sched_setaffinity(0,assigned))
        report['workers'].append(dict(host=host,gpu=gpu,device=device,cpu_affinity=assigned,
            controller_pid=proc.pid,process_start=Path(f'/proc/{proc.pid}/stat').read_text().rsplit(')',1)[1].split()[19],
            argv=argv,log=str(log),submitted=time.time()))
        atomic_write_json(receipt,report)
        time.sleep(.5)
    report['status']='controllers_submitted_actual_predictors_pending'
    atomic_write_json(receipt,report)
print('Submitted',len(suitable),'independently reserved',host,'workers for the nine-state private vinyl-monomer queue')
