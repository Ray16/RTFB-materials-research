"""Observe live predictor identity/resources and original gate termination locally."""
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
plan=json.loads((base/'uma_plan.json').read_text())
raw=subprocess.check_output(['nvidia-smi','--query-compute-apps=gpu_uuid,pid,used_memory','--format=csv,noheader,nounits'],text=True,timeout=20)
resident={int(r[1]):dict(uuid=r[0],memory_MiB=r[2]) for r in csv.reader(raw.splitlines(),skipinitialspace=True) if len(r)==3}
records=[];active_cpus=[];active_gpus=[];gates=[]
for path in sorted(base.glob(f'*_{host}_launch.json')):
    launch=json.loads(path.read_text())
    for worker in launch['workers']:
        stat=Path(f"/proc/{worker['gate_pid']}/stat")
        try:
            fields=stat.read_text().rsplit(')',1)[1].split()
            same=fields[19]==worker['process_start'];state=fields[0]
        except FileNotFoundError:same,state=False,None
        gates.append(dict(pid=worker['gate_pid'],same_start=same,state=state,terminal=not same or state=='Z'))
        for index in worker['tasks']:
            task=plan['tasks'][index];status=base/'tasks'/task['key']/'status.json'
            item=dict(key=task['key'],gpu=worker['gpu'],stage=launch['stage'])
            if not status.exists():
                item['status']='before_point'
            else:
                report=json.loads(status.read_text());item.update(status=report['status'],pid=report['pid'],error=report.get('error'))
                if report['status']=='computing':
                    process=psutil.Process(report['pid']);env=process.environ()
                    assert process.cpu_affinity()==worker['cpus']
                    assert env['CUDA_VISIBLE_DEVICES']==str(worker['gpu'])
                    assert all(env[k]=='1' for k in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS'))
                    if process.pid in resident:
                        assert resident[process.pid]['uuid']==worker['device']['uuid']
                        item['actual_predictor']=resident[process.pid]
                        active_cpus.extend(worker['cpus']);active_gpus.append(worker['gpu'])
            records.append(item)
assert len(active_cpus)==len(set(active_cpus)) and len(active_gpus)==len(set(active_gpus))
report=dict(status='actual_UMA_processes_observed',host=host,records=records,gates=gates,
    active_predictors=sum('actual_predictor' in r for r in records),
    completed=sum(r['status']=='gas_UMA_state_and_hydrogen_topology_verified' for r in records),
    failed=sum(r['status']=='failed_preserving_outputs' for r in records),load_1m=os.getloadavg()[0],checked=time.time())
atomic_write_json(base/f'{host}_observation.json',report)
print(json.dumps(report))
