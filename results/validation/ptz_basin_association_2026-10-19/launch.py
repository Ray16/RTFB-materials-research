"""Use every suitable free local GPU up to the two independent reciprocal inputs."""
import os,sys,json,time,subprocess
from pathlib import Path
from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows.resources import free_slots,gpu_suitability,host_capacity
base=Path(__file__).resolve().parent;assert json.loads((base/'preflight.json').read_text())['ok'];host=os.uname().nodename.split('.')[0];assert host=='lambda2'
plan=json.loads((base/'plan.json').read_text());tasks=[t for t in plan['tasks'] if not (base/('launch_'+t['target']+'.json')).exists()];slots=[]
for h,gpu in free_slots([host]):
 try:slots.append((gpu,gpu_suitability(gpu)))
 except RuntimeError:pass
assert slots or not tasks,'No suitable free GPU; no quantum work started'
host_capacity(2*min(len(tasks),len(slots)),minimum_ram_gib=64)
env=dict(os.environ,PYTHONPATH=str(ROOT/'src'));env['PATH']=str(Path(sys.executable).parent)+os.pathsep+env.get('PATH','')
for k in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS'):env[k]='2'
for (gpu,device),task in zip(slots,tasks):
 name=task['target'];out=base/(host+'_'+name);assert not out.exists();pre=host_capacity(2,minimum_ram_gib=64)
 argv=[sys.executable,'-u','-m','redox.workflows.gpu','--gpu',str(gpu),'--out',str(out),'--threads','2','--',sys.executable,'-u',str(base/'worker.py'),'--target',name]
 with (base/(name+'.log')).open('w') as log:p=subprocess.Popen(argv,env=env,stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,cwd=ROOT,start_new_session=True)
 record=dict(id='ptz_reciprocal_'+name,host=host,gpu=gpu,controller_pid=p.pid,process_start=Path(f'/proc/{p.pid}/stat').read_text().rsplit(')',1)[1].split()[19],argv=argv,resources=dict(pre,gpu=device),submitted=time.time(),scientific_outputs_verified=False)
 atomic_write_json(base/('launch_'+name+'.json'),record);print(name,p.pid,flush=True);time.sleep(1)
records=[json.loads(p.read_text()) for p in sorted(base.glob('launch_*.json'))]
atomic_write_json(base/'launch.json',dict(workers=records,pending=[t['target'] for t in plan['tasks'] if not (base/('launch_'+t['target']+'.json')).exists()]))
