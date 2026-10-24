"""Replace an unstarted native target after preserving/stopping its observer."""
from pathlib import Path
import fcntl,json,os,signal,socket,subprocess,sys,time
import psutil
from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api
BASE=ROOT/'results/validation/pmdi_occupancy_2026-10-23';history=BASE/'baseline_observer_history'
path=BASE/'baseline_coordinator.json';old=json.loads(path.read_text());assert old['host']==socket.gethostname().split('.')[0]
p=psutil.Process(old['pid']);assert Path(f'/proc/{p.pid}/stat').read_text().rsplit(')',1)[1].split()[19]==old['process_start'];assert 'coordinate_candidate_baseline.py' in ' '.join(p.cmdline())
(history/f'before_first_response_native_{time.time_ns()}.json').write_bytes(path.read_bytes())
p.send_signal(signal.SIGTERM)
for _ in range(100):
 if not p.is_running() or p.status()==psutil.STATUS_ZOMBIE:break
 time.sleep(.1)
else:raise RuntimeError('Observer did not stop; scientific workers remain untouched')
registry=BASE/'baseline_recoveries.json';r=json.loads(registry.read_text());key='baseline__design_pmdi_ethyl__occ_3';parent=Path(r['overrides'][key]['campaign'])
owner=api.load_module('first_response_owner_checks',ROOT/'scripts/pipeline/run_polymer_spacing.py')
for attempt in range(10):
 handles=[json.loads(q.read_text()) for q in parent.glob('thermal/launcher_handle*.json')]
 if not any(owner.owner_alive(h) for h in handles):break
 time.sleep(1)
else:raise RuntimeError('An earlier launcher remains live; keep its target until ownership resolves')
if any((parent/'thermal/claims').iterdir()) or any(parent.glob('thermal/workers/*.json')):
 raise RuntimeError('Earlier target acquired scientific work; preserve it instead of replacing')
new=ROOT/'results/polymer_models/campaigns/2026-10-23_crest_refinement/baseline__design_pmdi_ethyl__occ_3_native_first_response'
plan=api.read_plan(new);api.verify_sources(plan['producer_sha256']);api.verify_sources(plan['input_sha256']);records,invalid=api.admitted_records(new,'energy');assert len(records)==1 and not invalid
with registry.with_suffix('.lock').open('a+') as lock:
 fcntl.flock(lock,fcntl.LOCK_EX);r=json.loads(registry.read_text());assert Path(r['overrides'][key]['campaign'])==parent
 receipt=history/f'registry_before_first_response_native_{time.time_ns()}.json';receipt.write_bytes(registry.read_bytes())
 r['overrides'][key]=dict(campaign=str(new),hosts=['lambda2','lambda13'],reason='Exact-geometry fixed-density PCM response and both full spin controls passed; fresh native pilot with bounded first derivatives and contracted R; previous unstarted H100 fallback preserved')
 r['previous_registry']=str(receipt);r['updated']=time.time();atomic_write_json(registry,r)
remaining=(old['deadline']-time.time())/3600;assert 0<remaining<48
argv=[sys.executable,str(ROOT/'scripts/pipeline/coordinate_candidate_baseline.py'),'--campaign',str(BASE),'--hosts','lambda13','--max-hours',str(remaining),'--resume-observer']
log=history/f'first_response_resume_{time.time_ns()}.log'
with log.open('x') as f:
 child=subprocess.Popen(argv,cwd=ROOT,env=dict(os.environ,PYTHONPATH=str(ROOT/'src')),stdin=subprocess.DEVNULL,stdout=f,stderr=subprocess.STDOUT,start_new_session=True)
proof=dict(host=socket.gethostname(),pid=child.pid,process_start=Path(f'/proc/{child.pid}/stat').read_text().rsplit(')',1)[1].split()[19],argv=argv,started=time.time(),original_deadline=old['deadline'],log=str(log),earlier_observer_stopped=True,scientific_workers_preserved=True,old_unstarted_native_target=str(parent),new_native_target=str(new),registry_history=str(receipt))
atomic_write_json(history/'first_response_native_handoff.json',proof);print(json.dumps(proof))
