"""Read-only observation of the derivative and alternate-spin follow-up workers."""
import csv,json,os,subprocess,time
from pathlib import Path
import psutil
from redox.core.common import ROOT
from redox.core.io import atomic_write_json
cases=[('finer_derivative',ROOT/'results/validation/cata_cds_diagnosis_2026-10-19/finer_replay','cata_neutral/status.json'),('alternate_spin',ROOT/'results/validation/cata_spin_2026-10-19','parents/hyp_cata_ndi_ethyl/status.json')]
raw=subprocess.check_output(['nvidia-smi','--query-compute-apps=gpu_uuid,pid,used_memory','--format=csv,noheader,nounits'],text=True,timeout=20)
residents={int(r[1]):dict(uuid=r[0],used_memory_MiB=r[2]) for r in csv.reader(raw.splitlines(),skipinitialspace=True) if len(r)==3};reports=[]
for label,base,relative in cases:
 w=json.loads((base/'launch.json').read_text())['workers'][0];path=base/relative
 if not path.exists():reports.append(dict(case=label,status='worker_not_yet_reported'));continue
 r=json.loads(path.read_text());entry=dict(case=label,status=r['status'],pid=r['pid'],gpu=w['gpu'],error=r.get('error'),case_directory=str(base))
 try:
  p=psutil.Process(r['pid']);start=Path(f'/proc/{p.pid}/stat').read_text().rsplit(')',1)[1].split()[19];assert start==r['process_start']
  env=p.environ();aff=p.cpu_affinity();assert set(aff)==set(w['cpu_affinity'])
  assert env['CUDA_VISIBLE_DEVICES']==str(w['gpu']) and all(env[k]=='2' for k in ['OMP_NUM_THREADS','OMP_THREAD_LIMIT','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS'])
  entry.update(CPU_affinity=aff,native_thread_caps_verified=True,resident=residents.get(p.pid))
  if p.pid in residents:
   assert residents[p.pid]['uuid']==w['device']['uuid'];assert sum(v['uuid']==w['device']['uuid'] for v in residents.values())==1
   entry['actual_reserved_calculation_observed']=True
 except (psutil.NoSuchProcess,psutil.ZombieProcess):entry['process_live']=False
 reports.append(entry)
report=dict(status='followup_workers_observed',records=reports,load_1m=os.getloadavg()[0],checked=time.time(),physical_admission_not_inferred_from_startup=True)
atomic_write_json(ROOT/'results/validation/cata_solution_2026-10-19/followup_start_verification.json',report);print(json.dumps(report))
