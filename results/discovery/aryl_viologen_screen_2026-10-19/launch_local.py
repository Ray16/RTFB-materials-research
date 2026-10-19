"""Launch a disjoint declared subset on all suitable available local devices."""
import argparse
import csv
import json
import os
from pathlib import Path
import site
import socket
import subprocess
import sys
import time

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows.resources import GATE, allocation_guard, cpu_allocation, free_slots, host_capacity

base = Path(__file__).resolve().parent
parser = argparse.ArgumentParser()
parser.add_argument('--stage', choices=['pilot','screen'], required=True)
parser.add_argument('--tasks', help='disjoint comma-separated task indices for this host')
args = parser.parse_args()
host = socket.gethostname().split('.')[0]
plan = json.loads((base/'uma_plan.json').read_text())
indices = plan['pilot_indices'] if args.stage=='pilot' else list(map(int,args.tasks.split(',')))
assert len(indices)==len(set(indices)) and all(0<=i<len(plan['tasks']) for i in indices)
if args.stage=='screen':
    qualification=json.loads((base/'pilot_review.json').read_text())
    assert qualification['status']=='three_charge_70_atom_T4_resource_pilot_verified'
    assert not set(indices)&set(plan['pilot_indices'])
    assert qualification['plan_sha256']==__import__('hashlib').sha256((base/'uma_plan.json').read_bytes()).hexdigest()
for i in indices:
    assert not (base/'tasks'/plan['tasks'][i]['key']).exists(), 'Attempt already exists; inspect it'
preflight=json.loads((base/f'{host}_preflight.json').read_text())
assert preflight['ok'] and preflight['python']==sys.executable
receipt=base/f'{args.stage}_{host}_launch.json'
assert not receipt.exists(), 'Existing submission preserved; inspect original process identities'
(base/'claims').mkdir(exist_ok=True)
with allocation_guard():
    devices=[]
    for _, gpu in free_slots([host]):
        raw=subprocess.check_output(['nvidia-smi',f'--id={gpu}',
              '--query-gpu=index,uuid,name,memory.total,compute_cap','--format=csv,noheader,nounits'],text=True)
        index,uuid,name,memory,arch=next(csv.reader(raw.splitlines(),skipinitialspace=True))
        if args.stage=='pilot' and not (name=='Tesla T4' and arch=='7.5' and 14*1024<=float(memory)<=17*1024):
            continue
        if float(memory)<14*1024:
            continue
        devices.append(dict(index=int(index),uuid=uuid,name=name,memory_MiB=float(memory),architecture=arch))
    devices=devices[:len(indices)]
    assert devices
    resources=host_capacity(len(devices),minimum_ram_gib=12*len(devices))
    cpus=cpu_allocation(len(devices))
    env=dict(os.environ,PYTHONPATH=str(ROOT/'src'),FAIRCHEM_CACHE_DIR='/nfs/lambda_stor_01/homes/rzhu/.cache/fairchem')
    env['PATH']=str(Path(sys.executable).parent)+os.pathsep+env.get('PATH','')
    env['LD_LIBRARY_PATH']=':'.join(str(p) for b in site.getsitepackages() for p in Path(b).glob('nvidia/*/lib'))+':'+env.get('LD_LIBRARY_PATH','')
    for key in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
        env[key]='1'
    report=dict(status='submitting',stage=args.stage,host=host,resources=resources,workers=[],created=time.time())
    atomic_write_json(receipt,report)
    for i,device in enumerate(devices):
        assigned=indices[i::len(devices)]
        argv=[str(GATE),'run',str(device['index']),'--label','redox-aryl-UMA','--',
              'taskset','--cpu-list',str(cpus[i]),sys.executable,'-u',str(base/'worker.py'),
              '--tasks',','.join(map(str,assigned))]
        log=base/f"{args.stage}_{host}_gpu{device['index']}.log"
        with log.open('w') as stream:
            proc=subprocess.Popen(argv,cwd=ROOT,env=env,stdin=subprocess.DEVNULL,stdout=stream,
                stderr=subprocess.STDOUT,start_new_session=True,
                preexec_fn=lambda cpu=cpus[i]:os.sched_setaffinity(0,{cpu}))
        worker=dict(host=host,gpu=device['index'],device=device,cpus=[cpus[i]],tasks=assigned,
            gate_pid=proc.pid,process_start=Path(f'/proc/{proc.pid}/stat').read_text().rsplit(')',1)[1].split()[19],
            argv=argv,log=str(log),submitted=time.time())
        report['workers'].append(worker)
        atomic_write_json(receipt,report)
        time.sleep(.5)
    report['status']='controllers_submitted_actual_starts_pending'
    atomic_write_json(receipt,report)
print('Submitted',len(devices),'workers for',len(indices),'declared states on',host)
