"""Verify actual local queue owners, predictor descendants, CPUs and GPU identities."""
import csv
import json
import os
from pathlib import Path
import socket
import subprocess
import time

import psutil
from redox.core.io import atomic_write_json

base=Path(__file__).resolve().parent
host=socket.gethostname().split('.')[0]
path=base/f'launch_{host}.json'
assert path.exists(), 'Launch manifest not yet available; inspect original launch handle'
launch=json.loads(path.read_text())
raw=subprocess.check_output(['nvidia-smi','--query-compute-apps=gpu_uuid,pid,used_memory','--format=csv,noheader,nounits'],text=True,timeout=20)
resident={int(r[1]):dict(uuid=r[0],memory_MiB=r[2]) for r in csv.reader(raw.splitlines(),skipinitialspace=True) if len(r)==3}
records=[];cpus=[];gpus=[]
for worker in launch['workers']:
    status=base/'workers'/f"{host}_gpu{worker['gpu']}.json"
    item=dict(gpu=worker['gpu'],controller_pid=worker['controller_pid'],status='before_queue_owner')
    try:
        controller=psutil.Process(worker['controller_pid'])
        item['controller_process_state']=controller.status()
        item['controller_descendant_names']=[p.name() for p in controller.children(recursive=True)]
    except psutil.NoSuchProcess:
        item['controller_process_absent']=True
    if status.exists():
        owner=json.loads(status.read_text())
        item.update(queue_pid=owner['pid'],status=owner['status'],active_key=owner.get('active_key'),error=owner.get('error'))
        try:
            process=psutil.Process(owner['pid'])
            start=Path(f"/proc/{owner['pid']}/stat").read_text().rsplit(')',1)[1].split()[19]
            assert start==owner['process_start']
            env=process.environ()
            assert process.cpu_affinity()==worker['cpu_affinity']
            assert env['CUDA_VISIBLE_DEVICES']==str(worker['gpu'])
            assert all(env[k]=='2' for k in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS','XTB_THREADS'))
            tree=[process]+process.children(recursive=True)
            actual=[]
            for child in tree:
                if child.pid not in resident:continue
                assert resident[child.pid]['uuid']==worker['device']['uuid']
                assert child.cpu_affinity()==worker['cpu_affinity']
                actual.append(dict(pid=child.pid,argv=child.cmdline(),**resident[child.pid]))
            if actual:
                item.update(actual_GPU_processes=actual,cpu_affinity=process.cpu_affinity())
                cpus.extend(worker['cpu_affinity']);gpus.append(worker['gpu'])
            if item.get('active_key'):
                point=base/'tasks'/item['active_key']/'status.json'
                if point.exists():
                    state=json.loads(point.read_text());item['stage']=state.get('stage')
                    item['state_status']=state['status']
        except (psutil.NoSuchProcess,FileNotFoundError):
            item['original_queue_process_absent']=True
    execution=base/f"reserved_{host}_gpu{worker['gpu']}"/'execution.json'
    if execution.exists():
        r=json.loads(execution.read_text());item['controller_status']=r['status']
        item['controller_error']=r.get('error')
    records.append(item)
assert len(cpus)==len(set(cpus)) and len(gpus)==len(set(gpus))
states=[json.loads(p.read_text()) for p in (base/'tasks').glob('*/status.json')]
report=dict(status='actual_solution_queue_and_GPU_predictors_observed',host=host,workers=records,
    actual_GPU_workers=sum(bool(r.get('actual_GPU_processes')) for r in records),
    global_states_started=len(states),global_states_verified=sum(r['status']=='state_records_verified' for r in states),
    global_state_failures=sum(r['status']=='failed_preserved' for r in states),
    load_1m=os.getloadavg()[0],checked=time.time(),scientific_accuracy_validated=False)
atomic_write_json(base/f'observation_{host}.json',report)
print(json.dumps(report))
