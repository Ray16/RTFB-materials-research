"""Launch two fresh gas-reference cases on independent system reservations."""
import json,os,subprocess,sys,time
from pathlib import Path
from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows.resources import allocation_guard,cpu_allocation,free_slots,gpu_suitability,host_capacity
base=Path(__file__).resolve().parent
plan=json.loads((base/'plan.json').read_text());preflight=json.loads((base/'lambda13_preflight.json').read_text())
assert preflight['ok'] and len(plan['tasks'])==2
assert not (base/'launch.json').exists(),'Prior submission must be inspected, never overwritten'
with allocation_guard():
 resources=host_capacity(4,minimum_ram_gib=64);cpus=cpu_allocation(4)
 slots=free_slots(['lambda13']);assert len(slots)>=len(plan['tasks'])
 devices=[(host,gpu,gpu_suitability(gpu)) for host,gpu in slots[:len(plan['tasks'])]]
 env=dict(os.environ,PYTHONPATH=str(ROOT/'src'),FAIRCHEM_CACHE_DIR='/nfs/lambda_stor_01/homes/rzhu/.cache/fairchem')
 for key in ['OMP_NUM_THREADS','OMP_THREAD_LIMIT','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS','XTB_THREADS']:env[key]='2'
 record=dict(status='submitting',resources=resources,workers=[],source_plan=str(base/'plan.json'),actual_starts_require_verification=True)
 atomic_write_json(base/'launch.json',record)
 for i,(task,(host,gpu,device)) in enumerate(zip(plan['tasks'],devices)):
  assigned=cpus[i*2:i*2+2]
  argv=[sys.executable,'-u','-m','redox.workflows.gpu','--gpu',str(gpu),'--cpus',','.join(map(str,assigned)),'--threads','2','--out',str(base/f'reserved_{i}'),'--',sys.executable,'-u',str(base/'worker.py'),'--task',str(i)]
  with (base/f'controller_{i}.log').open('w') as output:
   child=subprocess.Popen(argv,cwd=ROOT,env=env,stdin=subprocess.DEVNULL,stdout=output,stderr=subprocess.STDOUT,start_new_session=True,preexec_fn=lambda assigned=assigned:os.sched_setaffinity(0,assigned))
  entry=dict(task=i,key=task['row']['id'],host=host,gpu=gpu,device=device,cpu_affinity=assigned,controller_pid=child.pid,process_start=Path(f'/proc/{child.pid}/stat').read_text().rsplit(')',1)[1].split()[19],argv=argv,submitted=time.time())
  record['workers'].append(entry);atomic_write_json(base/'launch.json',record)
 record['status']='two_controllers_submitted';atomic_write_json(base/'launch.json',record)
print('Submitted',len(record['workers']),'controllers')
