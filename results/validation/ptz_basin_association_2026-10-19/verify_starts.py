import json,subprocess,shlex,time
from pathlib import Path
root=Path('/nfs/lambda_stor_01/homes/rzhu/0_redox'); py='/nfs/lambda_stor_01/homes/rzhu/miniforge3/envs/redox/bin/python'
remote=r'''import json,sys,os,subprocess
from pathlib import Path
import psutil
jobs=json.loads(sys.argv[1]);out=[]
pools=('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')
r=subprocess.run(['nvidia-smi','--query-compute-apps=gpu_uuid,pid,used_memory','--format=csv,noheader,nounits'],capture_output=True,text=True,check=True)
resident={int(x.split(',')[1]):x.split(',')[0].strip() for x in r.stdout.splitlines() if x.strip()}
for j in jobs:
 p=psutil.Process(j['controller_pid']);start=Path(f'/proc/{p.pid}/stat').read_text().rsplit(')',1)[1].split()[19]
 assert start==j['process_start'],'PID identity changed'
 computed=[]
 for child in p.children(recursive=True):
  if child.pid not in resident:continue
  e=child.environ();aff=child.cpu_affinity();cmd=child.cmdline()
  assert len(aff)==2 and all(e.get(k)=='2' for k in pools)
  assert resident[child.pid]==j['resources']['gpu']['uuid']
  assert e.get('CUDA_VISIBLE_DEVICES')==str(j['gpu'])
  computed.append(dict(pid=child.pid,process_start=Path(f'/proc/{child.pid}/stat').read_text().rsplit(')',1)[1].split()[19],argv=cmd,cpu_affinity=aff,pools={k:e[k] for k in pools},gpu_uuid=resident[child.pid]))
 out.append(dict(id=j['id'],host=j['host'],gpu=j['gpu'],controller_pid=p.pid,resident_compute= computed,actual_GPU_start_verified=bool(computed)))
print(json.dumps(out))'''
byid={r['id']:r for r in json.loads((root/'results/validation/ptz_basin_association_2026-10-19/launch.json').read_text())['workers']}
results=[]
for host in sorted({x['host'] for x in byid.values()}):
 jobs=[j for j in byid.values() if j['host']==host]
 p=subprocess.run(['ssh','-o','BatchMode=yes',host,shlex.join([py,'-c',remote,json.dumps(jobs)])],capture_output=True,text=True,timeout=35)
 if p.returncode:raise RuntimeError(p.stderr)
 results.extend(json.loads(p.stdout))
target=root/'results/validation/ptz_basin_association_2026-10-19/start_verification.json';target.write_text(json.dumps(dict(created=time.time(),workers=results),indent=2)+'\n')
for r in results:print(r['host'],r['gpu'],r['id'],r['actual_GPU_start_verified'],[x['cpu_affinity'] for x in r['resident_compute']])
