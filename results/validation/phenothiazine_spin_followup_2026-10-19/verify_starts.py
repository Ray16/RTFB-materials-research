"""Verify actual reserved spin workers, rather than launcher or status files alone."""
import json
from pathlib import Path
import shlex
import subprocess
import time

base=Path(__file__).resolve().parent
py='/nfs/lambda_stor_01/homes/rzhu/miniforge3/envs/redox/bin/python'
program=r'''
import csv,json,os,psutil,subprocess,sys
from pathlib import Path
base=Path(sys.argv[1]);records=json.loads((base/'launch.json').read_text())['workers'];checked=[];used=set()
rows=list(csv.reader(subprocess.check_output(['nvidia-smi','--query-compute-apps=gpu_uuid,pid,used_memory','--format=csv,noheader,nounits'],text=True).splitlines(),skipinitialspace=True))
for r in records:
 control=psutil.Process(r['controller_pid'])
 assert Path(f'/proc/{control.pid}/stat').read_text().rsplit(')',1)[1].split()[19]==r['process_start']
 status=json.loads((base/'parents'/r['id']/'status.json').read_text());worker=psutil.Process(status['pid'])
 assert worker.pid in {p.pid for p in control.children(recursive=True)}
 assert Path(f'/proc/{worker.pid}/stat').read_text().rsplit(')',1)[1].split()[19]==status['process_start']
 expected=set(r['cpu_affinity']);assert expected==set(worker.cpu_affinity()) and not used&expected;used.update(expected)
 env=worker.environ()
 for name in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','XTB_THREADS'):assert env[name]=='2'
 assert env['CUDA_VISIBLE_DEVICES']==str(r['gpu'])
 resident=[dict(pid=int(p),memory_MiB=float(memory)) for uuid,p,memory in rows if uuid==r['device']['uuid']]
 assert len(resident)==1 and resident[0]['pid']==worker.pid
 assert status['status']=='computing'
 checked.append(dict(id=r['id'],gpu=r['gpu'],worker_pid=worker.pid,controller_pid=control.pid,
   cpu_affinity=sorted(expected),resident=resident,actual_quantum_start=True))
assert os.getloadavg()[0]<=os.cpu_count()
print(json.dumps(dict(host=os.uname().nodename.split('.')[0],workers=checked,load_1m=os.getloadavg()[0],online_cpus=os.cpu_count())))
'''
result=subprocess.run(['ssh','-o','BatchMode=yes','-o','ConnectTimeout=10','lambda1',
                       shlex.join([py,'-c',program,str(base)])],check=True,text=True,capture_output=True,timeout=45)
report=dict(status='actual_spin_starts_verified',observed=time.time(),**json.loads(result.stdout))
(base/'start_verification.json').write_text(json.dumps(report,indent=2)+'\n')
print(report['status'],len(report['workers']))
