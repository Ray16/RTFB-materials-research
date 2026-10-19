"""Start one reserved worker for each prepared state, using all needed free slots."""
import json
from pathlib import Path
import shlex
import subprocess
import time

root = Path('/nfs/lambda_stor_01/homes/rzhu/0_redox')
base = Path(__file__).resolve().parent
py = '/nfs/lambda_stor_01/homes/rzhu/miniforge3/envs/redox/bin/python'
plan = json.loads((base / 'plan.json').read_text())
assert len(plan['tasks']) == 12 and not (base / 'launch.json').exists()
for host in ('lambda1', 'lambda2'):
    preflight = json.loads((root / 'results/reviews/2026-10-19/ptz_benchmark_sources' /
                            (host + '_preflight.json')).read_text())
    assert preflight['ok'] and not preflight['errors']
program = r'''
import hashlib,json,os,subprocess,sys,time
from pathlib import Path
root=Path(sys.argv[1]);base=Path(sys.argv[2]);tasks=json.loads(sys.argv[3]);os.chdir(root)
sys.path.insert(0,str(root/'src'))
from redox.workflows.resources import allocation_guard,cpu_allocation,free_slots,gpu_suitability,host_capacity
from redox.qm.native_thermal import installed_versions
host=os.uname().nodename.split('.')[0]
assert installed_versions()=={'pyscf':'2.14.0','gpu4pyscf':'1.8.1'}
for path in (base/'plan.json',root/'library/manifest.csv',Path(sys.executable).parent/'xtb'):
 assert path.is_file(),str(path)
env=dict(os.environ,PYTHONPATH=str(root/'src'),REDOX_HESSIAN_MEMORY_GIB='6')
for key in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS','XTB_THREADS'):env[key]='2'
with allocation_guard():
 pre=host_capacity(2*len(tasks),minimum_ram_gib=8*len(tasks))
 cpus=cpu_allocation(2*len(tasks))
 free={g for h,g in free_slots([host])}
 for gpu,_ in tasks:assert gpu in free,'Selected GPU unavailable'
 devices={g:gpu_suitability(g) for g,_ in tasks}
 records=[]
 for offset,(gpu,task) in enumerate(tasks):
  selected=cpus[2*offset:2*offset+2]
  out=base/(host+'_task_'+str(task));assert not out.exists()
  argv=[sys.executable,'-u','-m','redox.workflows.gpu','--gpu',str(gpu),'--cpus',','.join(map(str,selected)),
        '--out',str(out),'--threads','2','--',sys.executable,'-u',str(base/'worker.py'),'--task',str(task)]
  log=base/(host+'_task_'+str(task)+'.log')
  with log.open('w') as output:
   child=subprocess.Popen(argv,cwd=root,env=env,stdout=output,stderr=subprocess.STDOUT,start_new_session=True)
  record=dict(host=host,gpu=gpu,task=task,cpu_affinity=selected,controller_pid=child.pid,
    process_start=Path(f'/proc/{child.pid}/stat').read_text().rsplit(')',1)[1].split()[19],
    argv=argv,log=str(log),output=str(out),device=devices[gpu],submitted=time.time(),
    actual_calculation_started=False,scientific_outputs_verified=False)
  records.append(record)
  (base/f'launch_{host}_{task}.json').write_text(json.dumps(record,indent=2)+'\n')
 # The child wrappers acquire this same guard after it is released. Fixed CPU
 # sets are disjoint; each wrapper rechecks them and atomically reserves its GPU.
(base/f'launch_{host}.json').write_text(json.dumps(dict(workers=records,host_preflight=pre),indent=2)+'\n')
print(host,len(records),'controllers submitted',flush=True)
'''
records = []
assignments = [('lambda1', [(i, i) for i in range(8)]),
               ('lambda2', [(i, i + 8) for i in range(4)])]
for host, tasks in assignments:
    argv = [py, '-c', program, str(root), str(base), json.dumps(tasks)]
    subprocess.run(['ssh', '-o', 'BatchMode=yes', '-o', 'ConnectTimeout=10', host,
                    shlex.join(argv)], check=True, timeout=60)
    records.extend(json.loads((base / f'launch_{host}.json').read_text())['workers'])
(base / 'launch.json').write_text(json.dumps(dict(status='controllers_submitted',
    workers=records, launched=time.time(), pending_tasks=[],
    actual_starts_require_verification=True, completed_scientific_states=0), indent=2) + '\n')
