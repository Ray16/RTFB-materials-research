"""Prepare a fully admitted cache; publish only under all molecule writer locks."""
import argparse,fcntl,hashlib,importlib.util,json,os,shutil,time
from contextlib import ExitStack
from pathlib import Path
from redox.core.common import ROOT,DFT
from redox.core.io import atomic_write_json
from redox.core.protocol import load_record,record_path
from redox.qm import native_thermal as N,hessian_memory as M
from redox.qm.projected_thermal import ensure_projected_thermal
from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity

parser=argparse.ArgumentParser();parser.add_argument('--adopt',action='store_true');args=parser.parse_args()
with allocation_guard():host_capacity(1);os.sched_setaffinity(0,cpu_allocation(1))
base=Path(__file__).resolve().parent;source=base/'finer_replay/cata_neutral/work';sd=DFT/'hyp_cata_adi_ethyl/neu';out=base/'admission';out.mkdir(exist_ok=True);prepared=out/'prepared_state'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
review_path=base/'finer_independent_review.json';review=json.loads(review_path.read_text());assert review['status']=='finer_target_independently_verified' and review['failed_checks']==0
for p,h in review['evidence_sha256'].items():assert sha(p)==h,p
spec=importlib.util.spec_from_file_location('independent_cache_admission',ROOT/'scripts/validation/audit/recompute_axes.py');audit=importlib.util.module_from_spec(spec);spec.loader.exec_module(audit)
def verify(directory):
 raw=json.loads((directory/'result.json').read_text());thermal=N.load_native_thermal(directory,0,1);assert thermal is not None
 correction,nimag,definition=audit.thermal_inputs(directory,raw)
 assert correction is not None and nimag==0 and definition==N.THERMAL_DEFINITION and abs(correction-thermal['g_thermal_total_eV'])<1e-8
 energy=load_record(directory,0,1);assert energy is not None
 return dict(native_thermal_eV=correction,native_G_eV=energy['e_smd_eV']+correction,minimum_frequency_cm=min(thermal['frequencies_cm']),maximum_asymmetry=thermal['maximum_asymmetry'])
receipt_path=out/'preparation.json'
existing_adoption=out/'adoption.json'
if existing_adoption.exists() and json.loads(existing_adoption.read_text()).get('status')=='reviewed_unchanged_geometry_cache_adopted':
 previous=json.loads(existing_adoption.read_text());receipt=json.loads(receipt_path.read_text());verify(sd)
 for name in ['opt.xyz','result.json']:assert sha(sd/name)==receipt['source_sha256'][name]
 native_path,_=N.cache_path(sd,0,1);assert sha(native_path)==receipt['prepared_sha256'][str(native_path.relative_to(sd))]
 assert all(sha(Path(previous['archive'])/p)==h for p,h in receipt['source_sha256'].items())
 print('Previously adopted cache and original archive verified; no write needed');raise SystemExit(0)
if not receipt_path.exists():
 old={str(p.relative_to(sd)):sha(p) for p in sd.rglob('*') if p.is_file()}
 assert sha(sd/'opt.xyz')==sha(source/'opt.xyz')
 shutil.copytree(sd,prepared)
 native,_=N.cache_path(source,0,1);target,_=N.cache_path(prepared,0,1)
 assert not target.exists(),'An existing addressed record needs a separate review'
 if target.with_suffix('').exists():shutil.rmtree(target.with_suffix(''))
 shutil.copytree(native.with_suffix(''),target.with_suffix(''))
 record=json.loads(native.read_text());qualification=json.loads((base/'finer_replay/cata_neutral/status.json').read_text())
 evidence={str(review_path):sha(review_path),**review['evidence_sha256']}
 record.update(status='derivative_checked_solution_thermal',adopted_from_isolated_target=str(source),admission_report=str(out/'adoption.json'))
 record['execution_derivatives']=dict(cds_step_Bohr=1e-6,qualification_only=False,production_admission=True,scientific_model_changed=False,qualification_sha256=M.derivative_certificate(N.NATIVE_THERMAL,1e-6),source_geometry_sha256=sha(sd/'opt.xyz'),worst_contact_indices=qualification['worst_contact_indices'],worst_contact_checks=qualification['worst_contact_checks'],target_evidence_sha256=evidence)
 atomic_write_json(target,record);N.select_policy(prepared)
 ensure_projected_thermal(target,prepared/'opt.xyz')
 values=verify(prepared)
 for name in ['opt.xyz','result.json']:assert sha(prepared/name)==sha(sd/name)
 assert sha(record_path(prepared,0,1))==sha(record_path(sd,0,1))
 assert old=={str(p.relative_to(sd)):sha(p) for p in sd.rglob('*') if p.is_file()}
 receipt=dict(status='qualified_cache_prepared_canonical_unchanged',canonical=str(sd),prepared=str(prepared),source_sha256=old,prepared_sha256={str(p.relative_to(prepared)):sha(p) for p in prepared.rglob('*') if p.is_file()},values=values,production_admission=False,created=time.time());atomic_write_json(receipt_path,receipt)
else:
 receipt=json.loads(receipt_path.read_text());verify(prepared)
print(json.dumps({k:receipt[k] for k in ['status','values','production_admission']}),flush=True)
if args.adopt:
 with ExitStack() as stack:
  try:
   for name in ['hyp_cata_adi_ethyl','hyp_cata_adi_ethyl_native_thermal','hyp_cata_adi_ethyl_spincheck']:
    lock=stack.enter_context((ROOT/'logs/pipeline_locks'/(name+'.lock')).open('a+'));fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
  except BlockingIOError:
   atomic_write_json(out/'adoption.json',dict(status='waiting_for_molecule_writers',canonical_modified=False,prepared_values=receipt['values'],updated=time.time()));print('Prepared cache retained; active molecule writers prevent adoption');raise SystemExit(0)
  assert receipt['source_sha256']=={str(p.relative_to(sd)):sha(p) for p in sd.rglob('*') if p.is_file()}
  assert receipt['prepared_sha256']=={str(p.relative_to(prepared)):sha(p) for p in prepared.rglob('*') if p.is_file()}
  archive=out/'superseded_canonical';assert not archive.exists();assert sd.stat().st_dev==out.stat().st_dev
  os.replace(sd,archive)
  try:
   os.replace(prepared,sd);values=verify(sd)
   assert all(sha(archive/p)==h for p,h in receipt['source_sha256'].items())
   for name in ['opt.xyz','result.json']:assert sha(sd/name)==sha(archive/name)
   assert sha(record_path(sd,0,1))==sha(record_path(archive,0,1))
  except BaseException:
   if sd.exists():os.replace(sd,out/'failed_install_preserved')
   os.replace(archive,sd);raise
  atomic_write_json(out/'adoption.json',dict(status='reviewed_unchanged_geometry_cache_adopted',canonical=str(sd),archive=str(archive),geometry_raw_and_active_energy_unchanged=True,values=values,source_sha256=receipt['source_sha256'],production_admission=True,completed=time.time()));print('Adopted under exclusive molecule/native/spin writer locks')
