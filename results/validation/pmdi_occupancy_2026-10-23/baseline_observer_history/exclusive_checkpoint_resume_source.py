from pathlib import Path
import json,os,subprocess,time,sys
import psutil
R=Path('/nfs/lambda_stor_01/homes/rzhu/0_redox');base=R/'results/validation/pmdi_occupancy_2026-10-23'
old=json.loads((base/'baseline_coordinator.json').read_text())
try:
 p=psutil.Process(old['pid']);start=Path(f'/proc/{p.pid}/stat').read_text().rsplit(')',1)[1].split()[19]
 if start==old['process_start']:raise RuntimeError('Observer is still live; do not duplicate it')
except psutil.NoSuchProcess:pass
remaining=(old['deadline']-time.time())/3600
assert 0<remaining<48
argv=[sys.executable,str(R/'scripts/pipeline/coordinate_candidate_baseline.py'),'--campaign',str(base),'--hosts','lambda13','--max-hours',str(remaining),'--resume-observer']
log=base/'baseline_observer_history'/f'exclusive_checkpoint_resume_{time.time_ns()}.log'
with log.open('x') as f:
 process=subprocess.Popen(argv,cwd=R,env=dict(os.environ,PYTHONPATH=str(R/'src')),stdin=subprocess.DEVNULL,stdout=f,stderr=subprocess.STDOUT,start_new_session=True)
start=Path(f'/proc/{process.pid}/stat').read_text().rsplit(')',1)[1].split()[19]
proof=dict(host='lambda2',pid=process.pid,process_start=start,argv=argv,log=str(log),started=time.time(),original_deadline=old['deadline'],scientific_workers_preserved=True)
(base/'baseline_observer_history/exclusive_checkpoint_observer_resume_launch.json').write_text(json.dumps(proof,indent=2)+'\n')
print(json.dumps(proof))
