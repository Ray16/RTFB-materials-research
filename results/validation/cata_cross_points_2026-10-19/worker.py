"""Compute disjoint, source-bound gas cross points; interpretation/QC is separate."""
import argparse,fcntl,hashlib,json,os,socket,time
from pathlib import Path
from redox.core.io import atomic_write_json
from redox.workflows.resources import verify_gate_owner
from redox.properties import reorg
base=Path(__file__).resolve().parent;parser=argparse.ArgumentParser();parser.add_argument('--tasks',required=True);a=parser.parse_args();plan=json.loads((base/'plan.json').read_text())
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
verify_gate_owner()
for p,h in plan['producer_sha256'].items():assert sha(p)==h,p
for index in map(int,a.tasks.split(',')):
 task=plan['tasks'][index];out=base/'tasks'/task['key'];out.mkdir(parents=True,exist_ok=True)
 lock=(out/'writer.lock').open('a+');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
 report=dict(status='starting',task=task,host=socket.gethostname(),pid=os.getpid(),process_start=Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()[19],CPU_affinity=sorted(os.sched_getaffinity(0)),started=time.time(),solution_kinetics_validated=False)
 def save():atomic_write_json(out/'status.json',report)
 save()
 try:
  for p,h in task['input_sha256'].items():assert sha(p)==h,p
  report.update(status='computing',physical_gpu=os.environ['CUDA_VISIBLE_DEVICES']);save()
  value=reorg._cross_energy_gas(task['id'],task['species'],task['charge'],task['mult'],task['at'],'gpu')
  path=reorg._cache_path(task['id'],task['species'],task['at'],task['charge'],task['mult']);record=json.loads(path.read_text())
  assert reorg._cache_ok(record,task['charge'],task['mult'],reorg._geom_hash(reorg.DFT/task['id']/task['at']/'opt.xyz'))
  for p,h in task['input_sha256'].items():assert sha(p)==h,p
  report.update(status='cross_point_solver_record_verified',energy_eV=value,record=str(path),record_sha256=sha(path),anion_unbound=record['anion_unbound'],gas_homo_eV=record['gas_homo_eV'])
 except BaseException as exc:
  report.update(status='failed_preserving_outputs',error=repr(exc));raise
 finally:
  report['finished']=time.time();save();lock.close()
