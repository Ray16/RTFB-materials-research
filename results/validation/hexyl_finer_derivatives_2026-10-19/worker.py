"""One source-frozen derivative replay and separate active-energy calculation."""
import argparse,hashlib,json,os,runpy,shutil,socket,sys,time,traceback
from pathlib import Path
from redox.core.io import atomic_write_json
from redox.core.protocol import load_record,record_path
from redox.core.common import smiles_to_xyz_map,read_manifest
from redox.qm.sp import ensure
from redox.workflows.resources import verify_gate_owner
base=Path(__file__).resolve().parent
a=argparse.ArgumentParser();a.add_argument('--task',required=True,type=int);arg=a.parse_args();verify_gate_owner()
plan=json.loads((base/'plan.json').read_text());task=plan['tasks'][arg.task]
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
for p,h in plan['producer_sha256'].items():assert sha(p)==h,p
status=dict(status='qualifying',task=task,host=socket.gethostname(),pid=os.getpid(),process_start=Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()[19],started=time.time(),production_admission=False,canonical_modified=False)
path=base/('worker_'+task['key']+'.json')
try:
 atomic_write_json(path,status)
 sys.argv=[str(base/'qualify.py'),'--manifest',str(base/'plan.json'),'--task',str(arg.task)]
 runpy.run_path(str(base/'qualify.py'),run_name='__main__')
 case=base/task['key'];qual=json.loads((case/'status.json').read_text());assert qual['status']=='isolated_target_verified'
 source=case/'work';energy_state=base/('energy_'+task['key']);energy_state.mkdir(exist_ok=False)
 for name in ['opt.xyz','result.json']:shutil.copy2(source/name,energy_state/name)
 row=next(r for r in read_manifest() if r['id']==task['id'] and r['state']==task['state'])
 assert smiles_to_xyz_map(row['smiles'],energy_state/'opt.xyz') is not None
 status.update(status='computing_active_energy',energy_state=str(energy_state));atomic_write_json(path,status)
 kind,sp=ensure(task['id'],task['state'],energy_state,task['charge'],task['mult'],backend='gpu');assert load_record(energy_state,task['charge'],task['mult'])
 assert all(sha(p)==h for p,h in qual['evidence_sha256'].items()),'Qualification evidence changed during SP'
 active=record_path(energy_state,task['charge'],task['mult'])
 status.update(status='private_native_and_SP_verified_pending_review',active_SP=str(active),active_SP_sha256=sha(active),e_smd_eV=sp['e_smd_eV'],minimum_frequency_cm=qual['minimum_frequency_cm'],finished=time.time())
except BaseException as exc:
 status.update(status='failed_preserved',error=repr(exc),traceback=traceback.format_exc(),finished=time.time());atomic_write_json(path,status);raise
atomic_write_json(path,status)
