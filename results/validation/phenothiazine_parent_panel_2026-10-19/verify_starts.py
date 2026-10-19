"""Inspect live owners, their actual GPU predictors and inherited CPU/thread caps."""
import json
from pathlib import Path
import shlex
import subprocess
import time

base = Path(__file__).resolve().parent
py = '/nfs/lambda_stor_01/homes/rzhu/miniforge3/envs/redox/bin/python'
launch = json.loads((base / 'launch.json').read_text())
program = r'''
import csv,json,os,subprocess,sys
from pathlib import Path
import psutil
records=json.loads(sys.argv[1]);report=[];used=set()
gpu_rows=list(csv.reader(subprocess.check_output(['nvidia-smi','--query-compute-apps=gpu_uuid,pid,used_memory','--format=csv,noheader,nounits'],text=True).splitlines(),skipinitialspace=True))
for record in records:
 execution=json.loads((Path(record['output'])/'execution.json').read_text())
 control=psutil.Process(record['controller_pid'])
 start=Path(f'/proc/{control.pid}/stat').read_text().rsplit(')',1)[1].split()[19]
 assert start==record['process_start']
 expected=set(record['cpu_affinity']);assert not expected&used;used.update(expected)
 descendants=control.children(recursive=True)
 predictors=[p for p in descendants if 'redox.qm.dft' in p.cmdline() or 'redox.qm.sp' in p.cmdline()]
 assert predictors,'No actual predictor under live controller'
 pidset={p.pid for p in descendants}
 resident=[dict(pid=int(row[1]),memory_MiB=float(row[2])) for row in gpu_rows if row[0]==record['device']['uuid']]
 assert resident and all(p['pid'] in pidset for p in resident),'Foreign process or missing GPU start'
 assert len(resident)==1,'Expected one resident quantum predictor per GPU'
 for process in predictors:
  assert set(process.cpu_affinity())==expected
  env=process.environ()
  for key in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS','XTB_THREADS'):assert env[key]=='2'
  assert env['CUDA_VISIBLE_DEVICES']==str(record['gpu'])
 report.append(dict(task=record['task'],gpu=record['gpu'],cpu_affinity=sorted(expected),
   controller_pid=control.pid,predictor_pids=[p.pid for p in predictors],resident=resident,
   actual_GPU_prediction_started=True,scientific_output_complete=False))
assert os.getloadavg()[0]<=os.cpu_count()
print(json.dumps(dict(host=os.uname().nodename.split('.')[0],load_1m=os.getloadavg()[0],online_cpus=os.cpu_count(),workers=report)))
'''
reports = []
for host in sorted({w['host'] for w in launch['workers']}):
    records = [w for w in launch['workers'] if w['host'] == host]
    argv = [py, '-c', program, json.dumps(records)]
    result = subprocess.run(['ssh', '-o', 'BatchMode=yes', '-o', 'ConnectTimeout=10', host,
                             shlex.join(argv)], check=True, text=True, capture_output=True, timeout=45)
    reports.append(json.loads(result.stdout))
receipt = dict(status='all_actual_starts_verified', observed=time.time(),
               actual_GPU_workers=sum(len(r['workers']) for r in reports),
               independent_scientific_admission=False, hosts=reports)
assert receipt['actual_GPU_workers'] == 12
(base / 'start_verification.json').write_text(json.dumps(receipt, indent=2) + '\n')
print(receipt['status'], receipt['actual_GPU_workers'])
