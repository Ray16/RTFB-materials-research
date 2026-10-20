"""Read-only proof of live reserved processes, GPU residency and disjoint CPUs."""
import json
from pathlib import Path
import shlex
import subprocess
import time

BASE=Path(__file__).resolve().parent
PYTHON='/nfs/lambda_stor_01/homes/rzhu/miniforge3/envs/redox/bin/python'
REMOTE=r'''
import csv,json,os,subprocess,time
from pathlib import Path
b=Path('/nfs/lambda_stor_01/homes/rzhu/0_redox/results/grafting/campaigns/2026-10-19_base_assisted_pilot')
host=os.uname().nodename.split('.')[0]
cards=list(csv.reader(subprocess.check_output(['nvidia-smi','--query-gpu=index,uuid','--format=csv,noheader,nounits'],text=True).splitlines(),skipinitialspace=True))
uuid={int(r[0]):r[1].strip() for r in cards}
resident=list(csv.reader(subprocess.check_output(['nvidia-smi','--query-compute-apps=gpu_uuid,pid,used_gpu_memory','--format=csv,noheader,nounits'],text=True).splitlines(),skipinitialspace=True))
result=[]
for p in sorted((b/'workers').glob(host+'_gpu*.json')):
 w=json.loads(p.read_text());pid=w['pid'];proc=Path('/proc')/str(pid)
 entry=dict(worker=w,live=False)
 try:
  start=(proc/'stat').read_text().rsplit(')',1)[1].split()[19]
  entry['live']=start==w['process_start']
  entry['actual_cpus']=sorted(os.sched_getaffinity(pid))
  env=dict(item.split('=',1) for item in (proc/'environ').read_bytes().decode().split('\0') if '=' in item)
  names=['OMP_NUM_THREADS','OMP_THREAD_LIMIT','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS','XTB_THREADS']
  entry['native_thread_limits']={n:env.get(n) for n in names}
  entry['visible_gpu']=env.get('CUDA_VISIBLE_DEVICES')
  parent=int((proc/'stat').read_text().rsplit(')',1)[1].split()[1])
  gate=(Path('/nfs/lambda_stor_01/homes/rzhu/.gpu_locks')/f'{host}_gpu{w["gpu"]}'/'owner').read_text().split()
  entry['gate_parent_matches']=gate[0]==str(parent) and gate[5]==(Path('/proc')/str(parent)/'stat').read_text().rsplit(')',1)[1].split()[19]
  entry['gpu_processes']=[dict(pid=int(r[1]),memory_MiB=int(r[2])) for r in resident if r[0].strip()==uuid[w['gpu']]]
  entry['own_gpu_resident']=any(r['pid']==pid for r in entry['gpu_processes'])
  entry['verified']=bool(entry['live'] and entry['actual_cpus']==w['cpus'] and all(env.get(n)=='2' for n in names) and entry['visible_gpu']==str(w['gpu']) and entry['gate_parent_matches'] and entry['own_gpu_resident'] and len(entry['gpu_processes'])==1)
 except (OSError,ValueError,IndexError) as exc:entry['error']=repr(exc);entry['verified']=False
 result.append(entry)
cpus=[set(r.get('actual_cpus',[])) for r in result]
overlap=any(a&c for i,a in enumerate(cpus) for c in cpus[i+1:])
print(json.dumps(dict(host=host,checked=time.time(),load_1m=os.getloadavg()[0],online_cpus=os.cpu_count(),CPU_overlap=overlap,workers=result)))
'''
reports=[]
for host in ('lambda13',):
 argv=['ssh','-n','-o','BatchMode=yes','-o','ConnectTimeout=8',host,shlex.join(['env','OMP_NUM_THREADS=1','OPENBLAS_NUM_THREADS=1',PYTHON,'-c',REMOTE])]
 proc=subprocess.run(argv,text=True,capture_output=True,timeout=55)
 if proc.returncode:reports.append(dict(host=host,error=proc.stderr,verified=False));continue
 reports.append(json.loads(proc.stdout))
out=dict(checked=time.time(),hosts=reports,verified_GPU_workers=sum(w.get('verified',False) for h in reports for w in h.get('workers',[])))
(BASE/'start_verification.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(dict(verified_GPU_workers=out['verified_GPU_workers'],hosts=[dict(host=h['host'],workers=len(h.get('workers',[])),overlap=h.get('CPU_overlap')) for h in reports])))
