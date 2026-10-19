"""Submit only the prerequisite restricted reproduction through the shared gate."""
import json,shlex,subprocess
from pathlib import Path
root=Path('/nfs/lambda_stor_01/homes/rzhu/0_redox');base=Path(__file__).resolve().parent
py='/nfs/lambda_stor_01/homes/rzhu/miniforge3/envs/redox/bin/python'
assert not (base/'launch_baseline.json').exists()
program=r'''import os,sys,json,subprocess,time
from pathlib import Path
root=Path(sys.argv[1]);base=Path(sys.argv[2]);os.chdir(root);sys.path.insert(0,str(root/'src'))
from redox.workflows.resources import free_slots,gpu_suitability,host_capacity
free=free_slots(['lambda13']);assert free
host_capacity(2,minimum_ram_gib=16);gpu=free[0][1];device=gpu_suitability(gpu)
env=dict(os.environ,PYTHONPATH=str(root/'src'))
for name in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS','XTB_THREADS'):env[name]='2'
argv=[sys.executable,'-u','-m','redox.workflows.gpu','--gpu',str(gpu),'--out',str(base/'reserved_baseline'),'--threads','2','--',sys.executable,'-u',str(base/'worker.py'),'--baseline']
with (base/'baseline_controller.log').open('w') as output:p=subprocess.Popen(argv,cwd=root,env=env,stdout=output,stderr=subprocess.STDOUT,stdin=subprocess.DEVNULL,start_new_session=True)
r=dict(host='lambda13',gpu=gpu,device=device,controller_pid=p.pid,process_start=Path(f'/proc/{p.pid}/stat').read_text().rsplit(')',1)[1].split()[19],argv=argv,submitted=time.time(),production_admission=False)
(base/'launch_baseline.json').write_text(json.dumps(r,indent=2)+'\n');print('baseline controller',p.pid,flush=True)
'''
subprocess.run(['ssh','-n','-o','BatchMode=yes','-o','ConnectTimeout=10','lambda13',shlex.join([py,'-c',program,str(root),str(base)])],check=True,timeout=45)
