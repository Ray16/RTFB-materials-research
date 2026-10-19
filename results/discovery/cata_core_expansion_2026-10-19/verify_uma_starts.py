"""Run on lambda13: verify actual reserved predictors, not submission alone."""
import csv,json,os,subprocess,time,hashlib
from pathlib import Path
import psutil
from redox.core.io import atomic_write_json
base=Path(__file__).resolve().parent
launch=json.loads((base/'uma_launch.json').read_text());records=[]
raw=subprocess.check_output(['nvidia-smi','--query-compute-apps=gpu_uuid,pid,used_memory','--format=csv,noheader,nounits'],text=True,timeout=20)
residents={int(row[1]):dict(uuid=row[0],used_memory_MiB=row[2]) for row in csv.reader(raw.splitlines(),skipinitialspace=True) if len(row)==3}
cpus=[];GPUs=[]
for task in launch['workers']:
 report=base/'uma_workers'/(task['id']+'__'+task['state'])/'status.json'
 if not report.exists():
  pid=task['controller_pid'];live=False
  try:live=Path(f'/proc/{pid}/stat').read_text().rsplit(')',1)[1].split()[19]==task['process_start']
  except FileNotFoundError:pass
  records.append(dict(task=task['task'],id=task['id'],state=task['state'],status='controller_alive_before_predictor' if live else 'controller_missing_before_predictor'));continue
 r=json.loads(report.read_text());entry=dict(task=task['task'],id=task['id'],state=task['state'],status=r['status'],pid=r['pid'])
 if r['status']=='gas_UMA_state_verified':
  assert hashlib.sha256(Path(r['result']).read_bytes()).hexdigest()==r['result_sha256']
  entry['output_sha256']=r['result_sha256'];entry['physical_gpu']=r['physical_gpu']
 elif r['status']=='computing':
  p=psutil.Process(r['pid']);start=Path(f'/proc/{p.pid}/stat').read_text().rsplit(')',1)[1].split()[19]
  env=p.environ();affinity=p.cpu_affinity()
  assert start==r['process_start'] and set(affinity)==set(task['cpus'])
  assert env['CUDA_VISIBLE_DEVICES']==str(task['gpu'])
  assert env['FAIRCHEM_CACHE_DIR']=='/nfs/lambda_stor_01/homes/rzhu/.cache/fairchem'
  assert all(env[k]=='2' for k in ['OMP_NUM_THREADS','OMP_THREAD_LIMIT','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS'])
  log=(base/f"uma_reserved_{task['task']}"/'command.log').read_text(errors='replace')
  entry.update(CPU_affinity=affinity,physical_gpu=env['CUDA_VISIBLE_DEVICES'],resident=residents.get(p.pid),actual_scan_started='[run ] '+task['id']+'/'+task['state'] in log,threads_verified=True)
  if p.pid in residents:
   assert residents[p.pid]['uuid']==task['device']['uuid']
   assert sum(v['uuid']==task['device']['uuid'] for v in residents.values())==1
  cpus.extend(affinity);GPUs.append(task['gpu'])
 records.append(entry)
assert len(cpus)==len(set(cpus)) and len(GPUs)==len(set(GPUs))
started=sum(bool(r['status']=='gas_UMA_state_verified' or (r.get('resident') and r.get('actual_scan_started'))) for r in records)
result=dict(status='six_actual_gas_UMA_starts_verified' if started==6 else 'partial_start_verification',actual_starts=int(started),completed=sum(r['status']=='gas_UMA_state_verified' for r in records),records=records,load_1m=os.getloadavg()[0],online_CPUs=os.cpu_count(),verified=time.time(),solution_potentials_completed=0,registered_candidates=0)
atomic_write_json(base/'uma_start_verification.json',result);print(json.dumps(result))
