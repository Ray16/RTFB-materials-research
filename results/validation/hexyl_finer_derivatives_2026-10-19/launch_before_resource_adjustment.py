"""Launch the two source-frozen independent branches on suitable free H100s."""
from pathlib import Path
import json,os,sys,subprocess,time
from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows.resources import free_slots,gpu_suitability,host_capacity
base=Path(__file__).resolve().parent;assert json.loads((base/'preflight.json').read_text())['ok']
assert os.uname().nodename.split('.')[0]=='lambda13';assert (base/'plan.json').is_file()
slots=free_slots(['lambda13']);gpus=[]
for host,index in slots:
 try:
  d=gpu_suitability(index)
  if 'H100' in d['name'] and d['memory_MiB']>=90000:gpus.append((index,d))
 except RuntimeError:continue
assert len(gpus)>=2,'Need two suitable free devices for two branches'
host_capacity(4,minimum_ram_gib=64)
env=dict(os.environ,PYTHONPATH=str(ROOT/'src'),REDOX_HESSIAN_MEMORY_GIB='24');env['PATH']=str(Path(sys.executable).parent)+os.pathsep+env.get('PATH','')
for k in ['OMP_NUM_THREADS','OMP_THREAD_LIMIT','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS']:env[k]='2'
workers=[]
for (gpu,device),branch in zip(gpus,['minus','plus']):
 key='hexyl_fine_'+branch;out=base/('lambda13_'+key);assert not out.exists();host_capacity(2,minimum_ram_gib=64)
 argv=[sys.executable,'-u','-m','redox.workflows.gpu','--gpu',str(gpu),'--out',str(out),'--threads','2','--',sys.executable,'-u',str(base/'worker.py'),'--task',str(['minus','plus'].index(branch))]
 with (base/(key+'.log')).open('w') as log:p=subprocess.Popen(argv,env=env,stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,cwd=ROOT,start_new_session=True)
 rec=dict(id=key,host='lambda13',gpu=gpu,controller_pid=p.pid,process_start=Path(f'/proc/{p.pid}/stat').read_text().rsplit(')',1)[1].split()[19],argv=argv,resources=dict(gpu=device),submitted=time.time(),scientific_outputs_verified=False)
 workers.append(rec);atomic_write_json(base/('launch_'+key+'.json'),rec);print(key,p.pid,flush=True);time.sleep(1)
atomic_write_json(base/'launch.json',dict(workers=workers))
