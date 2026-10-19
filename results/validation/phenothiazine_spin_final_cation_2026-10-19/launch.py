"""Launch ready molecule groups without reserving devices for pending sources."""
import json
from pathlib import Path
import shlex
import subprocess

base=Path(__file__).resolve().parent
root=base.parents[2]
py='/nfs/lambda_stor_01/homes/rzhu/miniforge3/envs/redox/bin/python'
assert not (base/'launch.json').exists()
preflight=json.loads((base/'lambda13_preflight.json').read_text())
assert preflight['ok'] and not preflight['errors']
program=r'''
import os,json,subprocess,sys,time,hashlib
from pathlib import Path
root=Path(sys.argv[1]);base=Path(sys.argv[2]);os.chdir(root);sys.path.insert(0,str(root/'src'))
from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity,free_slots,gpu_suitability
plan=json.loads((base/'plan.json').read_text());count=len(plan['parents']);records=[]
assert 0<count<=len({p['id'] for p in plan['parents']})
for p,h in plan['producer_sha256'].items():assert hashlib.sha256(Path(p).read_bytes()).hexdigest()==h,p
env=dict(os.environ,PYTHONPATH=str(root/'src'),REDOX_HESSIAN_MEMORY_GIB='6')
for key in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','XTB_THREADS'):env[key]='2'
host=os.uname().nodename.split('.')[0]
with allocation_guard():
 resources=host_capacity(2*count,minimum_ram_gib=8*count)
 cpus=cpu_allocation(2*count)
 free=free_slots([host]);devices=[]
 for _,gpu in free:
  try:devices.append(gpu_suitability(gpu))
  except RuntimeError:continue
 assert len(devices)>=count,'Not enough suitable free devices; no workers started'
 for index,device in enumerate(devices[:count]):
  gpu=device['index'];selected=cpus[2*index:2*index+2];out=base/(host+'_parent_'+str(index));assert not out.exists()
  argv=[sys.executable,'-u','-m','redox.workflows.gpu','--gpu',str(gpu),'--cpus',','.join(map(str,selected)),
        '--threads','2','--out',str(out),'--',sys.executable,'-u',str(base/'worker.py'),'--parent',str(index)]
  log=base/(host+'_parent_'+str(index)+'.log')
  with log.open('w') as stream:p=subprocess.Popen(argv,cwd=root,env=env,stdout=stream,stderr=subprocess.STDOUT,stdin=subprocess.DEVNULL,start_new_session=True)
  record=dict(host=host,gpu=gpu,parent=index,id=plan['parents'][index]['id'],device=device,
    cpu_affinity=selected,controller_pid=p.pid,process_start=Path(f'/proc/{p.pid}/stat').read_text().rsplit(')',1)[1].split()[19],
    argv=argv,output=str(out),log=str(log),submitted=time.time(),actual_GPU_start_verified=False)
  records.append(record);(base/f'launch_parent_{index}.json').write_text(json.dumps(record,indent=2)+'\n')
(base/'launch.json').write_text(json.dumps(dict(status='controllers_submitted',workers=records,resources=resources,
 current_ready_cases=sum(len(p['cases']) for p in plan['parents']),pending_primary_states=plan['pending_primary_states']),indent=2)+'\n')
print('controllers_submitted',len(records),flush=True)
'''
subprocess.run(['ssh','-n','-o','BatchMode=yes','-o','ConnectTimeout=10','lambda13',
                shlex.join([py,'-c',program,str(root),str(base)])],check=True,timeout=60)
