"""Bounded, source-addressed solution DFT alternate-spin validation.

The selected state's production SP must exist before a GPU is claimed. Helpers
are preserved; a lower alternative is promoted only after its native minimum,
connectivity and active energy are verified and the canonical writer is free.
"""
from __future__ import annotations
import fcntl,hashlib,json,math,os,shutil,time
from pathlib import Path
from redox.core.common import ROOT,DFT,UMA,read_manifest,smiles_to_xyz_map
from redox.core.protocol import ACTIVE_SP,load_record,record_path
from redox.qm.native_thermal import load_native_thermal,ensure_native_thermal,select_policy,cache_path,NATIVE_THERMAL
BASE=ROOT/'results/validation/physical_thermal_spin_2026-10-10/new_candidate_spin_checks'
NEAR_GAP_EV=.2


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def atomic_json(path,data):
 path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
 tmp=path.with_name(f'.{path.name}.{os.getpid()}.{time.time_ns()}.tmp')
 tmp.write_text(json.dumps(data,indent=2));tmp.replace(path)
def declared(gid,state):return next(r for r in read_manifest() if r['id']==gid and r['state']==state)
def snapshot(gid,state):
 sd=DFT/gid/state;p=sd/'result.json'
 if not p.exists() or not (sd/'opt.xyz').exists():return None
 raw=json.loads(p.read_text());q,m=int(raw['charge']),int(raw['mult']);sp=load_record(sd,q,m)
 if sp is None:return None
 row=declared(gid,state)
 return dict(id=gid,state=state,charge=q,mult=m,geometry_sha256=sha(sd/'opt.xyz'),
  smiles_sha256=hashlib.sha256(row['smiles'].encode()).hexdigest(),sp_record_sha256=sha(record_path(sd,q,m)),protocol=ACTIVE_SP)
def alternate_mult(m):
 if m in (1,3):return 4-m
 if m in (2,4):return 6-m
 raise ValueError('Only singlet/triplet or doublet/quartet comparisons are authorized here')
def comparison_path(gid,state,ident):
 key=hashlib.sha256(json.dumps(ident,sort_keys=True).encode()).hexdigest()[:16]
 return DFT/gid/state/'spin_validation'/f'{key}.json'
def read_comparison(gid,state):
 ident=snapshot(gid,state)
 if ident is None:return None
 p=comparison_path(gid,state,ident)
 if not p.exists():return None
 result=json.loads(p.read_text())
 if result.get('chosen_identity')!=ident:return None
 return result


def refreshed_electronic_gap(old_identity,new_identity,old_chosen,chosen,alternative,report):
 """Reuse an isolated spin comparison only after a like-for-like energy reproduction.

 Near-degenerate, thermochemical and promoted-spin cases require their full workflow.
 This function cannot turn a changed geometry, method or material energy change into
 a completed comparison. The caller must also prove the archived/recomputed records.
 """
 if {k:v for k,v in old_identity.items() if k!='sp_record_sha256'} != {
     k:v for k,v in new_identity.items() if k!='sp_record_sha256'}:
  raise ValueError('Spin comparison geometry, chemistry or protocol changed')
 if (report.get('status')!='spin_gap_verified' or report.get('thermochemical_gap_validated') is not False
     or report.get('ensemble_primary_integration_pending') or report.get('chosen_identity')!=old_identity):
  raise ValueError('Only completed isolated electronic comparisons can be refreshed')
 if (report.get('chosen_mult')!=new_identity['mult']
     or report.get('alternative_mult')!=alternative.get('mult')
     or alternative.get('charge')!=new_identity['charge']):
  raise ValueError('Spin comparison multiplicities or charge differ')
 for phase in ('smd','gas'):
  a,b=old_chosen.get(f'e_{phase}_eV'),chosen.get(f'e_{phase}_eV')
  if any(type(v) not in (int,float) or not math.isfinite(v) for v in (a,b)) or abs(a-b)>.001:
   raise ValueError('Chosen energy was not independently reproduced within 1 meV')
 values=(alternative.get('e_smd_eV'),alternative.get('smd_spin_contam'),report.get('electronic_gap_eV'))
 if any(type(v) not in (int,float) or not math.isfinite(v) for v in values) or abs(values[1])>.1:
  raise ValueError('Alternative spin diagnostics or old gap are invalid')
 if abs(values[0]-old_chosen['e_smd_eV']-values[2])>1e-8:
  raise ValueError('Original comparison no longer matches the helper energy')
 gap=values[0]-chosen['e_smd_eV']
 if min(gap,values[2])<=NEAR_GAP_EV:
  raise ValueError('Near-degenerate spin surfaces require matching thermochemistry')
 return gap


def refresh_reproduced_comparison(gid,state):
 """CPU-only derivation from a frozen independent SP reproduction and old report.

 All old evidence is retained. No quantum calculation or canonical state is changed.
 """
 from contextlib import ExitStack
 from redox.core.io import atomic_write_json
 from redox.workflows.revalidate_sp import done_valid,verified_candidate
 from redox.qm.dft import _opt_protocol_mismatch
 with ExitStack() as stack:
  for suffix in ('','_spincheck'):
   lock=stack.enter_context((ROOT/'logs/pipeline_locks'/f'{gid}{suffix}.lock').open('a+'))
   try:fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
   except BlockingIOError:return dict(id=gid,state=state,status='writer_busy')
  source=snapshot(gid,state)
  if source is None:return dict(id=gid,state=state,status='chosen_SP_missing')
  if read_comparison(gid,state):return dict(id=gid,state=state,status='current_comparison_present')
  sd=DFT/gid/state;path=record_path(sd,source['charge'],source['mult']);chosen=load_record(sd,source['charge'],source['mult'])
  evidence=chosen.get('independent_revalidation')
  if not evidence:return dict(id=gid,state=state,status='independent_reproduction_required')
  run=Path(evidence['run']);queue=json.loads((run/'queue.json').read_text())
  task=next(t for t in queue['tasks'] if t['key']==evidence['task'])
  reproduced=verified_candidate(task)
  if not done_valid(run,task) or reproduced is None or dict(chosen,independent_revalidation=None)!=dict(reproduced,independent_revalidation=None):
   raise ValueError('Current SP differs from the verified independent reproduction')
  if evidence['record_sha256']!=sha(record_path(Path(task['frozen']),source['charge'],source['mult'])):
   raise ValueError('Reproduction evidence hash changed')
  original=Path(task['frozen'])/'previous_record.json'
  if sha(original)!=evidence['original_sha256']:
   raise ValueError('Archived original SP hash changed')
  old_chosen=json.loads(original.read_text());matches=[]
  for old_path in (sd/'spin_validation').glob('*.json'):
   old_report=json.loads(old_path.read_text());old_source=old_report.get('chosen_identity') or {}
   if old_source.get('sp_record_sha256')==evidence['original_sha256']:matches.append((old_path,old_report,old_source))
  if len(matches)!=1:return dict(id=gid,state=state,status='unique_original_comparison_required')
  old_path,report,old_source=matches[0];helper=Path(report['helper']);raw=json.loads((helper/'result.json').read_text())
  if raw.get('spin_validation_status') not in (None,'validated') or _opt_protocol_mismatch(raw,source['charge']):
   raise ValueError('Alternative helper validation is incomplete')
  alt=load_record(helper,source['charge'],report['alternative_mult'])
  if alt is None:raise ValueError('Alternative active SP is invalid')
  if smiles_to_xyz_map(declared(gid,state)['smiles'],helper/'opt.xyz') is None:
   raise ValueError('Alternative helper changes molecular identity')
  gap=refreshed_electronic_gap(old_source,source,old_chosen,chosen,alt,report)
  old_sha=sha(old_path);alt_path=record_path(helper,source['charge'],report['alternative_mult']);alt_sha=sha(alt_path)
  derived=dict(report,chosen_identity=source,electronic_gap_eV=gap,created=time.time(),finished=time.time(),
               refreshed_from_path=str(old_path),refreshed_from_sha256=old_sha,
               independent_revalidation=evidence,alternative_record_path=str(alt_path),alternative_record_sha256=alt_sha,
               refresh_method='CPU recomputation of isolated electronic gap after independently reproduced chosen SP')
  if snapshot(gid,state)!=source or sha(old_path)!=old_sha or sha(alt_path)!=alt_sha:
   raise ValueError('Inputs changed while refreshing the comparison')
  target=comparison_path(gid,state,source);atomic_write_json(target,derived)
  return dict(id=gid,state=state,status='comparison_refreshed',path=str(target),electronic_gap_eV=gap)

def two_state_ensemble(chosen_G,alternative_G,temperature_K=298.15):
 if not all(isinstance(v,(int,float)) and math.isfinite(v) for v in (chosen_G,alternative_G,temperature_K)) or temperature_K<=0:
  raise ValueError('Finite energies and positive temperature required')
 kt=8.617333262145e-5*temperature_K;dg=alternative_G-chosen_G
 low=min(chosen_G,alternative_G);effective=low-kt*math.log1p(math.exp(-abs(dg)/kt))
 population=1/(1+math.exp(dg/kt)) if abs(dg/kt)<700 else float(dg<0)
 return dict(ensemble_free_energy_eV=effective,ensemble_correction_from_electronic_selected_eV=effective-chosen_G,alternative_population=population)


def _archive_tree(source,destination):
 shutil.copytree(source,destination)
 sums={str(p.relative_to(destination)):sha(p) for p in destination.rglob('*') if p.is_file()}
 assert all(sha(destination/name)==value for name,value in sums.items())
 atomic_json(destination.parent/(destination.name+'_SHA256.json'),sums)


def _helper_native(gid,state,sd):
 """Resolve a real helper saddle using at most two conservative signed-mode branches."""
 import numpy as np
 from ase.io import read,write
 from pyscf import gto
 from pyscf.hessian import thermo
 from pyscf.geomopt.geometric_solver import optimize
 from redox.qm.dft import _build_smd_mf,_ecp_for,_free_gpu_pool
 from redox.qm.sp import ensure,_ENERGY_KEYS
 try:return ensure_native_thermal(gid,state,state_dir=sd)[1]
 except RuntimeError as exc:
  if not ('not stationary' in str(exc) or 'not a minimum' in str(exc)):raise
 raw=json.loads((sd/'result.json').read_text());q,m=raw['charge'],raw['mult'];atoms=read(sd/'opt.xyz');symbols=atoms.get_chemical_symbols();n=len(atoms)
 folder=BASE/'helper_minimum_repairs'/gid/state/str(time.time_ns());folder.mkdir(parents=True)
 basis=NATIVE_THERMAL['anion_basis'] if q<0 else NATIVE_THERMAL['basis']
 def build(coords):
  mol=gto.M(atom=list(zip(symbols,coords.tolist())),basis=basis,ecp=_ecp_for(symbols,basis),charge=q,spin=m-1,verbose=0,max_memory=6000)
  mf=_build_smd_mf(mol,'r2scan','acetonitrile','gpu',disp='d4');mf.grids.level=3;mf.conv_tol=1e-12;mf.conv_tol_grad=1e-8;mf.max_cycle=200
  return mol,mf
 p,_=cache_path(sd,q,m);h=p.with_suffix('')/'hessian_cartesian.npy'
 if h.exists():
  mol,_=build(atoms.positions);H=np.load(h).reshape(n,3,n,3).transpose(0,2,1,3)
  modes=thermo.harmonic_analysis(mol,H,imaginary_freq=False);mode=np.asarray(modes['norm_mode'][int(np.argmin(modes['freq_wavenumber']))]);mode/=max(np.linalg.norm(mode,axis=1));seeds=[-.1*mode,.1*mode]
 else:seeds=[np.zeros_like(atoms.positions)]
 successes=[];evidence=[]
 for index,delta in enumerate(seeds):
  work=folder/f'branch_{index}';work.mkdir()
  try:
   _,mf=build(atoms.positions+delta);mf.kernel()
   if not mf.converged:raise RuntimeError('Helper refinement SCF failed')
   grad=mf.nuc_grad_method();grad.grid_response=True
   mol=optimize(grad.as_scanner(),maxsteps=120,assert_convergence=True,convergence_energy=1e-8,convergence_grms=3e-5,convergence_gmax=5e-5,convergence_drms=1e-4,convergence_dmax=2e-4)
   final=atoms.copy();final.positions=mol.atom_coords(unit='Angstrom');write(work/'opt.xyz',final);del mf,grad;_free_gpu_pool()
   if smiles_to_xyz_map(declared(gid,state)['smiles'],work/'opt.xyz') is None:raise RuntimeError('Helper refinement changes chemistry')
   r=dict(raw,optimized=True,g_thermal_eV=None)
   for key in _ENERGY_KEYS:r[key]=None
   atomic_json(work/'result.json',r);th=ensure_native_thermal(gid,state,state_dir=work)[1]
   _,energy=ensure(gid,state,work,q,m,backend='gpu')
   if load_record(work,q,m) is None:raise RuntimeError('Refined helper active SP failed')
   proto=energy['protocol'];r.update({key:energy.get(key) for key in _ENERGY_KEYS});r.update(sp_xc=proto['xc'],sp_basis=proto['basis'],sp_disp=proto['disp'],sp_nlc=proto['nlc'],xc=proto['xc'],basis=proto['basis'])
   atomic_json(work/'result.json',r);successes.append((energy['e_smd_eV'],work,th));evidence.append(dict(branch=index,status='verified',e_smd_eV=energy['e_smd_eV']))
  except Exception as exc:evidence.append(dict(branch=index,status='failed',error=repr(exc)))
  atomic_json(folder/'status.json',dict(branches=evidence,status='running'))
 if not successes:raise RuntimeError(f'Bounded helper minimum repair failed: {folder}')
 _,selected,th=min(successes,key=lambda x:x[0]);archive=folder/'superseded_helper';_archive_tree(sd,archive)
 for name in ['result.json','opt.xyz']:shutil.copy2(selected/name,sd/name)
 for name in ['sp','native_thermal']:
  (sd/name).mkdir(exist_ok=True)
  for source in (selected/name).iterdir():
   target=sd/name/source.name
   if target.exists():raise RuntimeError('Unexpected addressed helper collision; preserve both')
   if source.is_dir():shutil.copytree(source,target)
   else:shutil.copy2(source,target)
 select_policy(sd);atomic_json(folder/'status.json',dict(status='helper_minimum_verified',branches=evidence,selected=str(selected),archive=str(archive)))
 assert load_native_thermal(sd,q,m)
 return th


def _promote(gid,state,helper,source,report_path):
 """Returns False when the live motif writer still owns the canonical state."""
 lock=(ROOT/'logs/pipeline_locks'/f'{gid}.lock').open('a+')
 try:
  try:fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
  except BlockingIOError:return False
  if snapshot(gid,state)!=source:raise RuntimeError('Chosen source changed before spin promotion')
  old=DFT/gid/state;candidate=json.loads((helper/'result.json').read_text());q,m=candidate['charge'],candidate['mult']
  if not load_native_thermal(helper,q,m) or not load_record(helper,q,m):raise RuntimeError('Promotion requires validated own minimum and active SP')
  if smiles_to_xyz_map(declared(gid,state)['smiles'],helper/'opt.xyz') is None:raise RuntimeError('Promotion changes molecular connectivity')
  archive=BASE/'promotions'/gid/state/str(time.time_ns());archive.mkdir(parents=True);_archive_tree(old,archive/'canonical')
  reorg=DFT/gid/'reorg'
  if reorg.exists():_archive_tree(reorg,archive/'dependent_reorg')
  previous=ROOT/'calcs/spincheck'/gid/f'{state}_m{source["mult"]}'
  if previous.exists():_archive_tree(previous,archive/'previous_helper');shutil.rmtree(previous)
  shutil.copytree(old,previous)
  prepared=archive/'prepared_canonical';shutil.copytree(helper,prepared)
  r=json.loads((prepared/'result.json').read_text());r.update(id=gid,state=state,g_thermal_eV=None,n_imag=None,spin_promotion_archive=str(archive),spin_promotion_report=str(report_path));atomic_json(prepared/'result.json',r)
  select_policy(prepared)
  assert load_record(prepared,q,m) and load_native_thermal(prepared,q,m)
  # A complete prepared state is published by directory rename. Readers may briefly
  # see a missing state, but never a new geometry paired with old spin metadata.
  retired=archive/'retired_live_canonical';old.rename(retired)
  try:prepared.rename(old)
  except Exception:
   retired.rename(old);raise
  if reorg.exists():
   for path in reorg.glob('*.json'):path.unlink()
  override=ROOT/'config/spin_overrides.json'
  with (ROOT/'config/.spin_overrides.lock').open('a+') as ol:
   fcntl.flock(ol,fcntl.LOCK_EX);data=json.loads(override.read_text()) if override.exists() else {}
   data.setdefault(gid,{})[state]=dict(mult=m,charge=q,smiles_sha256=source['smiles_sha256'],protocol=ACTIVE_SP,
    validation='Verified lower active-level solution electronic spin with derivative-checked native minimum',source=str(report_path))
   atomic_json(override,data)
  return str(archive)
 finally:lock.close()


def check(gid,state):
 from redox.qm.dft import dft_smd,_opt_protocol_mismatch
 from redox.qm.sp import ensure
 source=snapshot(gid,state)
 if source is None:return dict(status='waiting_for_chosen_active_SP',id=gid,state=state)
 q,m=source['charge'],source['mult'];alt=alternate_mult(m);canonical=DFT/gid/state;chosen=load_record(canonical,q,m)
 helper=ROOT/'calcs/spincheck'/gid/f'{state}_m{alt}';run=BASE/'cases'/gid/state/str(time.time_ns());run.mkdir(parents=True)
 report=dict(id=gid,state=state,status='running',chosen_identity=source,chosen_mult=m,alternative_mult=alt,near_gap_threshold_eV=NEAR_GAP_EV,helper=str(helper),created=time.time(),host=os.uname().nodename,pid=os.getpid())
 path=run/'result.json';atomic_json(path,report)
 try:
  with (ROOT/'logs/pipeline_locks'/f'{gid}_spincheck.lock').open('a+') as lock:
   fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
   existing=json.loads((helper/'result.json').read_text()) if (helper/'result.json').exists() else None
   if helper.exists() and (existing is None or existing.get('spin_source_identity')!=source or _opt_protocol_mismatch(existing,q)):
    _archive_tree(helper,run/'superseded_helper');shutil.rmtree(helper);existing=None
   helper.mkdir(parents=True,exist_ok=True)
   if existing is None:
    result=dft_smd(canonical/'opt.xyz',q,alt,backend='gpu',opt_out=helper/'opt.xyz',do_freq=False)
    if result.get('converged_smd') is not True:raise RuntimeError('Alternative solution SCF did not converge')
    result.update(id=gid,state=f'{state}_m{alt}',parent_state=state,spin_source_identity=source,spin_validation_status='pending',g_thermal_eV=None)
    atomic_json(helper/'result.json',result)
   if smiles_to_xyz_map(declared(gid,state)['smiles'],helper/'opt.xyz') is None:raise RuntimeError('Alternative changes declared connectivity')
   _,energy=ensure(gid,f'{state}_m{alt}',helper,q,alt,backend='gpu');energy=load_record(helper,q,alt)
   if energy is None or abs(energy.get('smd_spin_contam',math.inf))>.1:raise RuntimeError('Alternative energy/spin QC failed')
   gap=energy['e_smd_eV']-chosen['e_smd_eV'];report.update(electronic_gap_eV=gap,thermochemical_gap_validated=False)
   if gap<-1e-4 or abs(gap)<=NEAR_GAP_EV:
    _helper_native(gid,state,helper);select_policy(helper);energy=load_record(helper,q,alt);gap=energy['e_smd_eV']-chosen['e_smd_eV'];report['electronic_gap_eV']=gap
   if snapshot(gid,state)!=source:
    raise RuntimeError('Chosen source changed during comparison; preserve alternative and requeue')
   if abs(gap)<=NEAR_GAP_EV:
    ct=load_native_thermal(canonical,q,m);at=load_native_thermal(helper,q,alt)
    if ct is None:
     report['status']='near_gap_waiting_chosen_native';atomic_json(path,report);atomic_json(comparison_path(gid,state,source),report);return report
    gc=chosen['e_smd_eV']+ct['g_thermal_total_eV'];ga=energy['e_smd_eV']+at['g_thermal_total_eV'];dg=ga-gc;kt=8.617333262145e-5*298.15
    ensemble=two_state_ensemble(gc,ga);effective=ensemble['ensemble_free_energy_eV'];palt=ensemble['alternative_population']
    report.update(thermochemical_gap_validated=True,thermochemical_gap_eV=dg,native_dependencies=dict(chosen=sha(cache_path(canonical,q,m)[0]),alternative=sha(cache_path(helper,q,alt)[0])),ensemble_free_energy_eV=effective,ensemble_correction_from_electronic_selected_eV=effective-gc,alternative_population=palt,single_surface_kinetics_validated=not (.01<palt<.99),ensemble_primary_integration_pending=True)
    assert ct['thermal_partition_protocol']==at['thermal_partition_protocol']
    report.update(native_partition_protocol=ct['thermal_partition_protocol'],native_partition_dependencies=dict(chosen=ct['thermal_postprocessing_sha256'],alternative=at['thermal_postprocessing_sha256']),native_partition_evidence_sha256={t['thermal_postprocessing_path']:t['thermal_postprocessing_sha256'] for t in (ct,at)})
    report['native_dependency_paths']=dict(chosen=str(cache_path(canonical,q,m)[0]),alternative=str(cache_path(helper,q,alt)[0]))
   if gap<-1e-4:
    archive=_promote(gid,state,helper,source,path)
    if archive is False:report['status']='lower_spin_waiting_motif_release';atomic_json(path,report);atomic_json(comparison_path(gid,state,source),report);return report
    report.update(status='lower_spin_promoted',archive=archive,old_chosen_identity=source,chosen_identity=snapshot(gid,state),chosen_mult=alt,alternative_mult=m,original_electronic_gap_eV=gap,electronic_gap_eV=-gap)
    previous=ROOT/'calcs/spincheck'/gid/f'{state}_m{m}'
    report.update(original_helper=str(helper),helper=str(previous))
    if report['thermochemical_gap_validated']:
     report.update(thermochemical_gap_eV=-report['thermochemical_gap_eV'],alternative_population=1-report['alternative_population'],ensemble_correction_from_electronic_selected_eV=effective-ga)
     report['native_dependencies']=dict(chosen=sha(cache_path(canonical,q,alt)[0]),alternative=sha(cache_path(previous,q,m)[0]))
     report['native_dependency_paths']=dict(chosen=str(cache_path(canonical,q,alt)[0]),alternative=str(cache_path(previous,q,m)[0]))
     current,prior=load_native_thermal(canonical,q,alt),load_native_thermal(previous,q,m)
     report.update(native_partition_dependencies=dict(chosen=current['thermal_postprocessing_sha256'],alternative=prior['thermal_postprocessing_sha256']),native_partition_evidence_sha256={t['thermal_postprocessing_path']:t['thermal_postprocessing_sha256'] for t in (current,prior)})
    previous_raw=json.loads((previous/'result.json').read_text())
    previous_raw.update(spin_validation_status='validated',spin_validation_report=str(path))
    atomic_json(previous/'result.json',previous_raw)
    from redox.properties.reorg import _couples,_cross_energy_gas
    for o,r in _couples(gid):
     if state not in (o[0],r[0]):continue
     report.setdefault('cross_refresh',[])
     for species,at in [(o,r[0]),(r,o[0])]:
      try:
       value=_cross_energy_gas(gid,species[0],species[1],species[2],at,'gpu')
       report['cross_refresh'].append(dict(species=species[0],at=at,status='computed',e_gas_eV=value))
      except Exception as exc:report['cross_refresh'].append(dict(species=species[0],at=at,status='failed',error=repr(exc)))
   else:report['status']='spin_gap_verified'
   raw=json.loads((helper/'result.json').read_text());raw['spin_validation_status']='validated';raw['spin_validation_report']=str(path);atomic_json(helper/'result.json',raw)
   report.update(finished=time.time(),experimental_accuracy_validated=False)
   atomic_json(path,report);atomic_json(comparison_path(gid,state,report['chosen_identity']),report);return report
 except Exception as exc:
  report.update(status='source_stale' if 'source changed' in str(exc).lower() else 'failed',error=repr(exc),finished=time.time());atomic_json(path,report)
  if (helper/'result.json').exists():
   raw=json.loads((helper/'result.json').read_text());raw['spin_validation_status']=report['status'];atomic_json(helper/'result.json',raw)
  return report
