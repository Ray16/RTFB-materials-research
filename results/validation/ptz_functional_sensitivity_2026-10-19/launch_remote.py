"""Use all free suitable GPUs for the nine fixed-geometry method points."""
import json,os,subprocess,sys,time
from pathlib import Path
from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows.resources import allocation_guard,cpu_allocation,free_slots,gpu_suitability,host_capacity
base=Path(__file__).resolve().parent;plan=json.loads((base/'plan.json').read_text());assert not (base/'launch.json').exists()
assert json.loads((base/'lambda13_preflight.json').read_text())['ok']
with allocation_guard():
 slots=[]
 for host,gpu in free_slots(['lambda13']):
  try:slots.append((host,gpu,gpu_suitability(gpu)))
  except RuntimeError:pass
 slots=slots[:len(plan['tasks'])];assert slots
 resources=host_capacity(2*len(slots),minimum_ram_gib=16*len(slots));cpus=cpu_allocation(2*len(slots))
 env=dict(os.environ,PYTHONPATH=str(ROOT/'src'))
 for k in ['OMP_NUM_THREADS','OMP_THREAD_LIMIT','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS','XTB_THREADS']:env[k]='2'
 record=dict(status='submitting',resources=resources,workers=[],points=len(plan['tasks']),scientific_QC_separate=True,production_admission=False)
 atomic_write_json(base/'launch.json',record)
 for i,(host,gpu,device) in enumerate(slots):
  assigned=cpus[i*2:i*2+2];indices=list(range(i,len(plan['tasks']),len(slots)))
  argv=[sys.executable,'-u','-m','redox.workflows.gpu','--gpu',str(gpu),'--cpus',','.join(map(str,assigned)),'--threads','2','--out',str(base/f'reserved_{i}'),'--',sys.executable,'-u',str(base/'worker.py'),'--tasks',','.join(map(str,indices))]
  with (base/f'controller_{i}.log').open('w') as log:child=subprocess.Popen(argv,cwd=ROOT,env=env,stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,start_new_session=True,preexec_fn=lambda assigned=assigned:os.sched_setaffinity(0,assigned))
  record['workers'].append(dict(host=host,gpu=gpu,device=device,tasks=indices,cpu_affinity=assigned,controller_pid=child.pid,process_start=Path(f'/proc/{child.pid}/stat').read_text().rsplit(')',1)[1].split()[19],argv=argv,submitted=time.time()));atomic_write_json(base/'launch.json',record)
 record['status']='all_available_suitable_method_point_controllers_submitted';atomic_write_json(base/'launch.json',record)
print('Submitted',len(record['workers']),'workers for',len(plan['tasks']),'points')
