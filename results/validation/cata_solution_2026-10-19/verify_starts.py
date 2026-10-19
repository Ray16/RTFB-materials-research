"""Observe actual molecular predictors and inherited allocations on lambda13."""
import csv,json,os,subprocess,time
from pathlib import Path
import psutil
from redox.core.io import atomic_write_json
base=Path(__file__).resolve().parent;launch=json.loads((base/'launch.json').read_text())
raw=subprocess.check_output(['nvidia-smi','--query-compute-apps=gpu_uuid,pid,used_memory','--format=csv,noheader,nounits'],text=True,timeout=20)
residents={int(r[1]):dict(uuid=r[0],used_memory_MiB=r[2]) for r in csv.reader(raw.splitlines(),skipinitialspace=True) if len(r)==3}
records=[];observed_cpus=[];observed_GPUs=[]
for worker in launch['workers']:
 p=base/'tasks'/worker['key']/'status.json';entry=dict(key=worker['key'],gpu=worker['gpu'],controller_pid=worker['controller_pid'])
 if not p.exists():
  stat=Path(f"/proc/{worker['controller_pid']}/stat")
  entry['status']='live_controller_before_worker' if stat.exists() and stat.read_text().rsplit(')',1)[1].split()[19]==worker['process_start'] else 'controller_missing'
 else:
  r=json.loads(p.read_text());entry.update(status=r['status'],stage=r.get('stage'),pid=r['pid'],error=r.get('error'))
  if r['status']=='state_records_verified':entry['completed_records']=True
  elif r['status'] in ('preparing','computing'):
   process=psutil.Process(r['pid']);start=Path(f'/proc/{process.pid}/stat').read_text().rsplit(')',1)[1].split()[19]
   assert start==r['process_start']
   tree=[process]+process.children(recursive=True);predictors=[]
   for child in tree:
    try:
     env=child.environ();affinity=child.cpu_affinity()
     assert set(affinity)==set(worker['cpu_affinity'])
     assert env.get('CUDA_VISIBLE_DEVICES')==str(worker['gpu'])
     assert all(env.get(k)=='2' for k in ['OMP_NUM_THREADS','OMP_THREAD_LIMIT','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS','XTB_THREADS'])
     if child.pid in residents:
      assert residents[child.pid]['uuid']==worker['device']['uuid']
      predictors.append(dict(pid=child.pid,argv=child.cmdline(),resident=residents[child.pid],CPU_affinity=affinity))
    except (psutil.NoSuchProcess,psutil.ZombieProcess):continue
   entry['actual_predictors']=predictors
   if predictors:
    assert sum(v['uuid']==worker['device']['uuid'] for v in residents.values())==len(predictors)==1
    observed_cpus.extend(worker['cpu_affinity']);observed_GPUs.append(worker['gpu'])
 assert entry.get('status')!='controller_missing',entry
 records.append(entry)
assert len(observed_cpus)==len(set(observed_cpus)) and len(observed_GPUs)==len(set(observed_GPUs))
started=sum(bool(r.get('actual_predictors') or r.get('completed_records')) for r in records)
report=dict(status='six_actual_solution_starts_verified' if started==6 else 'partial_solution_start_verification',actual_starts=started,complete_state_records=sum(bool(r.get('completed_records')) for r in records),records=records,host='lambda13',load_1m=os.getloadavg()[0],online_CPUs=os.cpu_count(),verified=time.time(),physical_admission_is_separate=True)
atomic_write_json(base/'start_verification.json',report);print(json.dumps(report))
