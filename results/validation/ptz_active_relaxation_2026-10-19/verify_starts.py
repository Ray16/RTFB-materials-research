"""Read-only resident-process and actual optimizer-progress verification."""
import csv,json,os,subprocess,time
from pathlib import Path
import psutil
from redox.core.io import atomic_write_json
base=Path(__file__).resolve().parent;plan=json.loads((base/'plan.json').read_text());launch=json.loads((base/'launch.json').read_text())
raw=subprocess.check_output(['nvidia-smi','--query-compute-apps=gpu_uuid,pid,used_memory','--format=csv,noheader,nounits'],text=True,timeout=20);resident={int(r[1]):dict(uuid=r[0],memory_MiB=r[2]) for r in csv.reader(raw.splitlines(),skipinitialspace=True) if len(r)==3};records=[];cpus=[];gpus=[]
for w in launch['workers']:
 state=plan['tasks'][w['task']]['state'];path=base/state/'status.json';entry=dict(state=state,gpu=w['gpu'],controller_pid=w['controller_pid'])
 if not path.exists():entry['status']='controller_before_predictor'
 else:
  r=json.loads(path.read_text());entry.update(status=r['status'],pid=r['pid'],cycle=r.get('cycle'),error=r.get('error'))
  if r['status']=='active_level_stationary_geometry_verified':entry['completed']=True
  elif r['status'] not in ('failed_preserved',):
   p=psutil.Process(r['pid']);start=Path(f'/proc/{p.pid}/stat').read_text().rsplit(')',1)[1].split()[19];assert start==r['process_start']
   env=p.environ();aff=p.cpu_affinity();assert set(aff)==set(w['cpu_affinity'])
   assert env['CUDA_VISIBLE_DEVICES']==str(w['gpu']) and all(env[k]=='2' for k in ['OMP_NUM_THREADS','OMP_THREAD_LIMIT','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS'])
   if p.pid in resident:
    assert resident[p.pid]['uuid']==w['device']['uuid'];entry.update(actual_predictor=resident[p.pid],CPU_affinity=aff,actual_optimizer_progress=bool(r.get('cycle')));cpus.extend(aff);gpus.append(w['gpu'])
 records.append(entry)
assert len(cpus)==len(set(cpus)) and len(gpus)==len(set(gpus))
started=sum(bool(r.get('completed') or (r.get('actual_predictor') and r.get('actual_optimizer_progress'))) for r in records)
report=dict(status='three_actual_private_optimization_starts_verified' if started==3 else 'partial_private_start_verification',actual_starts=started,records=records,load_1m=os.getloadavg()[0],online_CPUs=os.cpu_count(),verified=time.time(),production_admission=False)
atomic_write_json(base/'start_verification.json',report);print(json.dumps(report))
