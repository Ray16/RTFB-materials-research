#!/usr/bin/env python
"""Single reserved, resumable native Hessian worker for declared affected ladders."""
import fcntl,hashlib,json,os,socket,time
from pathlib import Path
from redox.core.common import ROOT,DFT,read_manifest,read_result
from redox.qm.native_thermal import ensure_native_thermal,select_policy,NATIVE_THERMAL,cache_path
from redox.properties.reorg import _couples,_cross_energy_gas
BASE=ROOT/'results/validation/physical_thermal_spin_2026-10-10'
STATUS=BASE/'native_physical_cohort_status.json'
def main():
 plan=json.loads((BASE/'native_physical_cohort_manifest.json').read_text())
 status=dict(status='running',host=socket.gethostname(),pid=os.getpid(),process_start=Path('/proc/self/stat').read_text().split()[21],
   cpus=sorted(os.sched_getaffinity(0)),gpu=os.getenv('CUDA_VISIBLE_DEVICES'),started=time.time(),items=[],cross_items=[])
 assert status['cpus']==[0,1] and status['gpu'] and os.getloadavg()[0]+2<os.cpu_count()
 def save():
  tmp=STATUS.with_suffix('.tmp');tmp.write_text(json.dumps(status,indent=2));tmp.replace(STATUS)
 save();failures=0
 for key in plan['states']:
  gid,state=key.split(':');sd=DFT/gid/state
  if not(sd/'result.json').exists():continue
  select_policy(sd)
  item=dict(id=gid,state=state,status='running',started=time.time(),source_sha256=hashlib.sha256((sd/'opt.xyz').read_bytes()).hexdigest())
  status['items'].append(item);save();print('NATIVE_START',key,flush=True)
  try:
   with (ROOT/'logs/pipeline_locks'/f'{gid}_native_thermal.lock').open('a+') as lock:
    fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    try:kind,record=ensure_native_thermal(gid,state)
    except RuntimeError as exc:
     repairable='not stationary' in str(exc) or 'not a minimum' in str(exc)
     if 'Partial native thermal outputs' in str(exc):
      for p in (BASE/'native_minimum_repairs'/gid/state).glob('*/status.json'):
       evidence=json.loads(p.read_text())
       if evidence.get('source_sha256')==item['source_sha256'] and ('not stationary' in evidence.get('reason','') or 'not a minimum' in evidence.get('reason','')):
        repairable=True
     if not repairable:raise
     from redox.qm.native_thermal_repair import repair_minimum
     kind,record=repair_minimum(gid,state,str(exc))
   item.update(status='verified',kind=kind,g_thermal_eV=record['g_thermal_total_eV'],record=str(cache_path(sd,record['identity']['charge'],record['identity']['mult'])[0]),finished=time.time());failures=0
   print('NATIVE_VERIFIED',key,record['g_thermal_total_eV'],flush=True)
  except Exception as exc:
   item.update(status='failed',error=repr(exc),finished=time.time());failures+=1;print('NATIVE_FAILED',key,repr(exc),flush=True)
  save()
  if failures>=3:status['status']='circuit_breaker';save();return
 # Restore only cross caches affected by previously verified canonical changes.
 for gid in ['ferrocene','d3_80izcq']:
  for o,r in _couples(gid):
   for (species,q,m),at in [(o,r[0]),(r,o[0])]:
    try:
     value=_cross_energy_gas(gid,species,q,m,at,'gpu');status['cross_items'].append(dict(id=gid,species=species,at=at,e_gas_eV=value,status='computed_or_present'))
    except Exception as exc:status['cross_items'].append(dict(id=gid,species=species,at=at,status='failed',error=repr(exc)))
    save()
 status.update(status='completed_with_failures' if any(i['status']=='failed' for i in status['items']+status['cross_items']) else 'completed',finished=time.time());save()
if __name__=='__main__':main()
