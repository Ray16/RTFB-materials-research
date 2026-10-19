"""Inspect actual cross-point predictors, completed records and inherited resources."""
import csv,json,os,subprocess,time
from pathlib import Path
import psutil
from redox.core.io import atomic_write_json
base=Path(__file__).resolve().parent;plan=json.loads((base/'plan.json').read_text());launch=json.loads((base/'launch.json').read_text())
raw=subprocess.check_output(['nvidia-smi','--query-compute-apps=gpu_uuid,pid,used_memory','--format=csv,noheader,nounits'],text=True,timeout=20);resident={int(r[1]):dict(uuid=r[0],memory_MiB=r[2]) for r in csv.reader(raw.splitlines(),skipinitialspace=True) if len(r)==3}
records=[];cpus=[];gpus=[]
for w in launch['workers']:
 for index in w['tasks']:
  task=plan['tasks'][index];p=base/'tasks'/task['key']/'status.json';entry=dict(key=task['key'],gpu=w['gpu'])
  if not p.exists():entry['status']='controller_before_case'
  else:
   r=json.loads(p.read_text());entry.update(status=r['status'],pid=r['pid'],error=r.get('error'))
   if r['status']=='fixed_geometry_electronic_method_point_verified':entry['completed_record']=str(p)
   elif r['status']=='computing':
    process=psutil.Process(r['pid']);env=process.environ();aff=process.cpu_affinity();assert set(aff)==set(w['cpu_affinity'])
    assert env['CUDA_VISIBLE_DEVICES']==str(w['gpu']) and all(env[k]=='2' for k in ['OMP_NUM_THREADS','OMP_THREAD_LIMIT','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS'])
    if process.pid in resident:
     assert resident[process.pid]['uuid']==w['device']['uuid'];entry.update(actual_predictor=resident[process.pid],CPU_affinity=aff)
     cpus.extend(aff);gpus.append(w['gpu'])
  records.append(entry)
assert len(cpus)==len(set(cpus)) and len(gpus)==len(set(gpus))
report=dict(status='functional_method_progress_observed',records=records,actual_started_or_complete=sum(bool(r.get('actual_predictor') or r.get('completed_record')) for r in records),completed=sum(bool(r.get('completed_record')) for r in records),active_predictors=sum(bool(r.get('actual_predictor')) for r in records),host='lambda13',load_1m=os.getloadavg()[0],verified=time.time(),solution_kinetics_not_validated=True)
atomic_write_json(base/'start_verification.json',report);print(json.dumps(report))
