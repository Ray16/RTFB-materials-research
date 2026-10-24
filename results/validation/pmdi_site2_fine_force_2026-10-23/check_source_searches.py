"""Read-only process/evidence check for the independent source searches."""
from pathlib import Path
import json,os,socket,time
import psutil
R=Path('/nfs/lambda_stor_01/homes/rzhu/0_redox')
host=socket.gethostname().split('.')[0]
records=[]
roots=[R/'calcs/polymer_models/2026-10-23_conformer_revision',R/'calcs/polymer_models/pmdi_occupancy_2026-10-23']
for root in roots:
 for path in sorted(root.glob('*/status.json')):
  state=json.loads(path.read_text())
  if state.get('status')!='running' or state.get('host')!=host:continue
  parent,predictor=psutil.Process(state['pid']),psutil.Process(state['predictor_pid'])
  for proc,field in ((parent,'process_start'),(predictor,'predictor_start')):
   assert Path(f'/proc/{proc.pid}/stat').read_text().rsplit(')',1)[1].split()[19]==state[field]
  assert predictor.ppid()==parent.pid
  assert predictor.cpu_affinity()==state['CPU_affinity']
  assert 'crest' in predictor.name().lower()
  artifacts={}
  for name in ('crest.log','crest_ensemble.xyz','ensemble_energies.log'):
   p=path.parent/name
   if p.exists():artifacts[name]=dict(bytes=p.stat().st_size,mtime=p.stat().st_mtime)
  records.append(dict(case=path.parent.name,worker_pid=parent.pid,predictor_pid=predictor.pid,predictor_process_start=state['predictor_start'],CPU_affinity=predictor.cpu_affinity(),predictor_CPU_time_seconds=sum(predictor.cpu_times()[:2]),artifacts=artifacts,search_status='running_not_final',DFT_reranking_complete=False))
proof=dict(checked=time.time(),host=host,verified_live_predictors=len(records),records=records)
out=R/'results/validation/pmdi_site2_fine_force_2026-10-23'/f'live_searches_{host}.json';out.write_text(json.dumps(proof,indent=2)+'\n')
print(json.dumps(dict(host=host,verified_live_predictors=len(records),completed_searches_in_this_record=0)))
