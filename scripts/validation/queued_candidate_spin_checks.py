#!/usr/bin/env python
"""CPU coordinator claims a GPU only for ready, source-addressed spin checks."""
import argparse,fcntl,hashlib,json,os,signal,socket,subprocess,sys,time
from pathlib import Path
from redox.core.common import ROOT,UMA,read_manifest,DFT
from redox.validation.spins import snapshot,check,read_comparison,atomic_json,BASE
GATE='/nfs/lambda_stor_01/homes/rzhu/bin/gpu_reserve'
SKILL='/nfs/lambda_stor_01/homes/rzhu/.codex-personal/skills/lambda-fleet/scripts/preflight_node.py'

def needed(ids,explicit=None):
 rows=read_manifest();tasks=[]
 for row in rows:
  gid,state=row['id'],row['state'];key=gid+':'+state
  if explicit is not None and key not in explicit and gid not in explicit:continue
  if explicit is None and gid not in ids:continue
  p=UMA/gid/state/'result.json'
  if not p.exists():continue
  uma=json.loads(p.read_text());mult=uma.get('chosen_mult');gap=uma.get('spin_gap_eV');scanned=uma.get('mults_scanned',[])
  even=mult in (1,3)
  eligible=(even and (gap is None or gap<.5 or mult!=int(row.get('mult_hint') or 1))) or (mult in (2,4) and (mult==4 or (4 in scanned and (gap is None or gap<.5))))
  if eligible:tasks.append(key)
 return tasks

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--worker',action='store_true');ap.add_argument('--only');ap.add_argument('--fleet',type=Path,default=ROOT/'results/validation/candidate_fleet_today_2026-10-11')
 ap.add_argument('--gpu',type=int,default=1);ap.add_argument('--cpus',default='2,3');ap.add_argument('--max-hours',type=float,default=12);a=ap.parse_args()
 BASE.mkdir(parents=True,exist_ok=True);explicit=set(a.only.split(',')) if a.only else None
 if a.worker:
  assert os.getenv('CUDA_VISIBLE_DEVICES') and sorted(os.sched_getaffinity(0))==[int(x) for x in a.cpus.split(',')]
  if explicit is None:ap.error('worker requires exact --only')
  results=[]
  for key in sorted(explicit):
   gid,state=key.split(':');result=check(gid,state);results.append(result);print('CANDIDATE_SPIN_RESULT',key,result['status'],result.get('electronic_gap_eV'),flush=True)
  atomic_json(BASE/'last_worker_result.json',dict(results=results,host=socket.gethostname(),pid=os.getpid(),process_start=Path('/proc/self/stat').read_text().split()[21],finished=time.time()))
  if any(r['status'] in ('failed','source_stale') for r in results):raise SystemExit(1)
  return
 status=dict(status='waiting_for_ready_chosen_SP',host=socket.gethostname(),pid=os.getpid(),process_start=Path('/proc/self/stat').read_text().split()[21],cpus=sorted(os.sched_getaffinity(0)),gpu=None,compute_gpu=a.gpu,compute_cpus=a.cpus,started=time.time(),max_hours=a.max_hours,attempts={},items=[])
 path=BASE/'coordinator_status.json';deadline=time.time()+a.max_hours*3600;child=None
 def save():atomic_json(path,status)
 def stop(sig,frame):raise KeyboardInterrupt(str(sig))
 signal.signal(signal.SIGTERM,stop);save()
 env=dict(os.environ,PYTHONPATH=str(ROOT/'src'));env['PATH']=str(Path(sys.executable).parent)+os.pathsep+env.get('PATH','')
 import site
 env['LD_LIBRARY_PATH']=':'.join(str(p) for base in site.getsitepackages() for p in Path(base).glob('nvidia/*/lib'))+':'+env.get('LD_LIBRARY_PATH','')
 try:
  while time.time()<deadline:
   queue=json.loads((a.fleet/'queue.json').read_text());tasks=needed(set(queue['ids']),explicit);status['discovered_tasks']=tasks;ready=[]
   for key in tasks:
    gid,state=key.split(':');source=snapshot(gid,state)
    if source is None:continue
    result=read_comparison(gid,state)
    if result and result['status'] in ('spin_gap_verified','lower_spin_promoted'):continue
    if result and result['status']=='near_gap_waiting_chosen_native':
     from redox.qm.native_thermal import load_native_thermal
     if load_native_thermal(DFT/gid/state,source['charge'],source['mult']) is None:continue
    if result and result['status']=='lower_spin_waiting_motif_release':
     with (ROOT/'logs/pipeline_locks'/f'{gid}.lock').open('a+') as guard:
      try:fcntl.flock(guard,fcntl.LOCK_EX|fcntl.LOCK_NB)
      except BlockingIOError:continue
    key_source=key+'__'+source['geometry_sha256']
    if status['attempts'].get(key_source,0)>=2:continue
    ready.append((key,key_source,source))
   status['ready_count']=len(ready);save()
   if not ready:
    status['status']='waiting_for_ready_chosen_SP_or_native';save();time.sleep(30);continue
   free=subprocess.check_output([GATE,'list','--hosts',socket.gethostname()],text=True).splitlines()
   if socket.gethostname()+' '+str(a.gpu) not in free:status['status']='waiting_for_shared_GPU';save();time.sleep(30);continue
   key,key_source,source=ready[0]
   claim=BASE/'claims'/key.replace(':','__');claim.parent.mkdir(exist_ok=True)
   try:claim.mkdir()
   except FileExistsError:raise RuntimeError('Existing spin claim ownership must be audited, not reclaimed by age')
   atomic_json(claim/'owner.json',dict(host=socket.gethostname(),pid=os.getpid(),process_start=status['process_start'],source=source))
   try:
    # Fresh dependency/CPU preflight every exact item; no GPU reservation held while pending.
    preflight=BASE/'preflights'/f'{key.replace(":","__")}_{time.time_ns()}.json';preflight.parent.mkdir(exist_ok=True)
    argv=['taskset','-c',a.cpus,sys.executable,SKILL,'local','--python',sys.executable,'--module','pyscf','--module','gpu4pyscf','--module','ase','--module','rdkit','--cpu-set',a.cpus,'--path',str(ROOT),'--min-available-ram-mib','10000']
    with preflight.open('w') as output:pre=subprocess.run(argv,stdout=output,stderr=subprocess.STDOUT)
    if pre.returncode:raise RuntimeError('Fresh spin dependency/CPU preflight failed')
    log=BASE/'logs'/f'{key.replace(":","__")}_{time.time_ns()}.log';log.parent.mkdir(exist_ok=True)
    argv=[GATE,'run',str(a.gpu),'--label','redox-source-addressed-candidate-spin-check','--','taskset','-c',a.cpus,'timeout','--signal=TERM','--kill-after=30','3h',sys.executable,'-u',str(Path(__file__).resolve()),'--worker','--only',key,'--cpus',a.cpus]
    entry=dict(key=key,source=source,argv=argv,preflight=str(preflight),log=str(log),started=time.time(),producer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    status['items'].append(entry);status['status']='running_actual_ready_spin_check';save()
    for name in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS'):env[name]='2'
    with log.open('w') as output:
     child=subprocess.Popen(argv,env=env,cwd=ROOT,stdout=output,stderr=subprocess.STDOUT,start_new_session=True);entry['gate_pid']=child.pid;save();rc=child.wait();child=None
    entry.update(returncode=rc,finished=time.time());status['attempts'][key_source]=status['attempts'].get(key_source,0)+1
    status['consecutive_failures']=status.get('consecutive_failures',0)+1 if rc else 0;save()
    if status['consecutive_failures']>=3:raise RuntimeError('Three consecutive actual spin-check failures; circuit breaker')
   finally:shutil_rmtree_claim(claim)
  status.update(status='deadline_reached',finished=time.time());save()
 except (Exception,KeyboardInterrupt) as exc:status.update(status='failed' if not isinstance(exc,KeyboardInterrupt) else 'terminated',error=repr(exc),finished=time.time());save();raise
 finally:
  if child is not None and child.poll() is None:
   os.killpg(child.pid,signal.SIGTERM)
   try:child.wait(timeout=30)
   except subprocess.TimeoutExpired:os.killpg(child.pid,signal.SIGKILL);child.wait()

def shutil_rmtree_claim(claim):
 import shutil
 shutil.rmtree(claim)
if __name__=='__main__':main()
