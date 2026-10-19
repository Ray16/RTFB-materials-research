"""Freeze two terminal coarse-CDS failures for the qualified 1e-6-Bohr replay."""
from pathlib import Path
import fcntl,hashlib,json,os,shutil,time
from redox.core.common import ROOT,read_manifest,smiles_to_xyz_map
from redox.core.io import atomic_write_json
from redox.qm.native_thermal import NATIVE_THERMAL
b=Path(__file__).resolve().parent;gid='exp_ndi_hexyl_pentyloxy';state='red1'
old=ROOT/'results/validation/physical_thermal_spin_2026-10-10/native_minimum_repairs'/gid/state/'1792377821041388905'
status=json.loads((old/'status.json').read_text());assert status['status']=='repair_failed' and len(status['branches'])==2
assert all('Hessian symmetry' in v['error'] for v in status['branches'])
for suffix in ('','_native_thermal'):
 with (ROOT/'logs/pipeline_locks'/f'{gid}{suffix}.lock').open('a+') as lock:fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
qualifier=b/'qualify.py';shutil.copy2(ROOT/'results/validation/ndi_methyl_radical_minimum_2026-10-19/qualify.py',qualifier)
row=next(r for r in read_manifest() if r['id']==gid and r['state']==state);tasks=[]
for label,sign in [('minus',-1),('plus',1)]:
 source=old/f'branch_{sign}';dest=b/('frozen_'+label);dest.mkdir(exist_ok=False)
 for name in ['opt.xyz','result.json']:shutil.copy2(source/name,dest/name)
 raw=json.loads((dest/'result.json').read_text());assert raw['charge']==-1 and raw['mult']==2
 assert smiles_to_xyz_map(row['smiles'],dest/'opt.xyz') is not None
 tasks.append(dict(key=label,id=gid,state=state,source=str(dest),charge=-1,mult=2,control=False,
  input_sha256={str(dest/n):sha(dest/n) for n in ['opt.xyz','result.json']},
  original_inputs_sha256={str(source/n):sha(source/n) for n in ['opt.xyz','result.json']},
  source_optimizer_convergence_not_claimed=True))
producers={str(qualifier):sha(qualifier)}
for relative in ['qm/native_thermal.py','qm/dft.py','qm/sp.py','qm/hessian_memory.py','core/common.py','core/protocol.py','qm/qrrho.py','qm/xtb_kernel.py','qm/projected_thermal.py']:
 p=ROOT/'src/redox'/relative;producers[str(p)]=sha(p)
plan=dict(protocol=NATIVE_THERMAL,producer_sha256=producers,tasks=tasks,cds_step_Bohr=1e-6,
 physical_thresholds_unchanged=True,memory_cap_GiB=24,maximum_attempts_per_input=1,
 canonical_modified=False,old_failure_status_path=str(old/'status.json'),old_failure_status_sha256=sha(old/'status.json'),created=time.time())
atomic_write_json(b/'plan.json',plan);print('Prepared two changed-derivative qualifications; no canonical writes')
