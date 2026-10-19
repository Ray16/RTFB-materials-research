"""Three independent frozen-state active-gradient controls on reserved H100s."""
import json,subprocess,shlex
from pathlib import Path
root=Path('/nfs/lambda_stor_01/homes/rzhu/0_redox');base=Path(__file__).resolve().parent
py='/nfs/lambda_stor_01/homes/rzhu/miniforge3/envs/redox/bin/python'
assert not (base/'launch.json').exists()
assert json.loads((base/'lambda13_preflight.json').read_text())['ok']
program=r'''import os,sys,json,hashlib,subprocess,time
from pathlib import Path
root=Path(sys.argv[1]);base=Path(sys.argv[2]);os.chdir(root);sys.path.insert(0,str(root/'src'))
from redox.workflows.resources import allocation_guard,cpu_allocation,free_slots,gpu_suitability,host_capacity
plan=json.loads((base/'plan.json').read_text());count=len(plan['tasks']);workers=[]
for path,h in plan['producer_sha256'].items():assert hashlib.sha256(Path(path).read_bytes()).hexdigest()==h,path
env=dict(os.environ,PYTHONPATH=str(root/'src'))
for key in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','XTB_THREADS'):env[key]='2'
with allocation_guard():
 resources=host_capacity(2*count,minimum_ram_gib=16*count);cpus=cpu_allocation(2*count)
 free=free_slots(['lambda13']);assert len(free)>=count,'No controls started: insufficient free devices'
 devices=[gpu_suitability(g) for _,g in free[:count]]
 for index,device in enumerate(devices):
  selected=cpus[2*index:2*index+2];gpu=device['index'];out=base/f'reserved_{index}';assert not out.exists()
  argv=[sys.executable,'-u','-m','redox.workflows.gpu','--gpu',str(gpu),'--cpus',','.join(map(str,selected)),
    '--threads','2','--out',str(out),'--',sys.executable,'-u',str(base/'worker.py'),'--task',str(index)]
  with (base/f'controller_{index}.log').open('w') as log:p=subprocess.Popen(argv,cwd=root,env=env,stdout=log,stderr=subprocess.STDOUT,stdin=subprocess.DEVNULL,start_new_session=True)
  record=dict(host='lambda13',task=index,state=plan['tasks'][index]['state'],gpu=gpu,device=device,CPU_affinity=selected,
    controller_pid=p.pid,process_start=Path(f'/proc/{p.pid}/stat').read_text().rsplit(')',1)[1].split()[19],argv=argv,submitted=time.time())
  workers.append(record);(base/f'launch_{index}.json').write_text(json.dumps(record,indent=2)+'\n')
(base/'launch.json').write_text(json.dumps(dict(status='controllers_submitted',workers=workers,resources=resources,production_admission=False),indent=2)+'\n')
print('gradient controllers',len(workers),flush=True)
'''
subprocess.run(['ssh','-n','-o','BatchMode=yes','-o','ConnectTimeout=10','lambda13',shlex.join([py,'-c',program,str(root),str(base)])],check=True,timeout=45)
