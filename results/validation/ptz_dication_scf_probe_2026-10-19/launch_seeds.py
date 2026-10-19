"""Run the four independent orbital-seed probes on separate suitable free H100s."""
import json,shlex,subprocess
from pathlib import Path
root=Path('/nfs/lambda_stor_01/homes/rzhu/0_redox');base=Path(__file__).resolve().parent
py='/nfs/lambda_stor_01/homes/rzhu/miniforge3/envs/redox/bin/python'
assert not (base/'launch_seeds.json').exists()
program=r'''import os,sys,json,hashlib,subprocess,time
from pathlib import Path
root=Path(sys.argv[1]);base=Path(sys.argv[2]);os.chdir(root);sys.path.insert(0,str(root/'src'))
from redox.workflows.resources import allocation_guard,cpu_allocation,free_slots,gpu_suitability,host_capacity
control=json.loads((base/'baseline/report.json').read_text());assert control['status']=='restricted_baseline_reproduced'
assert hashlib.sha256((base/'baseline/orbitals.npz').read_bytes()).hexdigest()==control['orbitals_sha256']
plan=json.loads((base/'plan.json').read_text());count=len(plan['seed_angles_radians']);records=[]
env=dict(os.environ,PYTHONPATH=str(root/'src'))
for name in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','XTB_THREADS'):env[name]='2'
with allocation_guard():
 resources=host_capacity(2*count,minimum_ram_gib=16*count);cpus=cpu_allocation(2*count)
 slots=free_slots(['lambda13']);assert len(slots)>=count,'Insufficient free devices; no probe started'
 devices=[gpu_suitability(g) for _,g in slots[:count]]
 for index,device in enumerate(devices):
  selected=cpus[2*index:2*index+2];gpu=device['index'];out=base/f'reserved_seed_{index}';assert not out.exists()
  argv=[sys.executable,'-u','-m','redox.workflows.gpu','--gpu',str(gpu),'--cpus',','.join(map(str,selected)),
    '--threads','2','--out',str(out),'--',sys.executable,'-u',str(base/'worker.py'),'--seed',str(index)]
  with (base/f'seed_controller_{index}.log').open('w') as log:p=subprocess.Popen(argv,cwd=root,env=env,stdout=log,stderr=subprocess.STDOUT,stdin=subprocess.DEVNULL,start_new_session=True)
  r=dict(host='lambda13',seed=index,gpu=gpu,device=device,CPU_affinity=selected,controller_pid=p.pid,
    process_start=Path(f'/proc/{p.pid}/stat').read_text().rsplit(')',1)[1].split()[19],argv=argv,submitted=time.time())
  records.append(r);(base/f'launch_seed_{index}.json').write_text(json.dumps(r,indent=2)+'\n')
receipt=dict(status='seed_controllers_submitted',workers=records,resources=resources,
 baseline_report_sha256=hashlib.sha256((base/'baseline/report.json').read_bytes()).hexdigest(),orbitals_sha256=control['orbitals_sha256'],production_admission=False)
(base/'launch_seeds.json').write_text(json.dumps(receipt,indent=2)+'\n');print('seed controllers',len(records),flush=True)
'''
subprocess.run(['ssh','-n','-o','BatchMode=yes','-o','ConnectTimeout=10','lambda13',shlex.join([py,'-c',program,str(root),str(base)])],check=True,timeout=45)
