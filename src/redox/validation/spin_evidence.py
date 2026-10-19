"""Read-only spin publication gates shared by candidate and scorecard readers.

These predicates do not run calculations or alter producer records.
"""

def spin_comparison_running(gid, lock_root=None):
 """Observe the actual writer lock, including while an old helper is archived."""
 import fcntl
 from redox.core.common import ROOT
 directory=ROOT/'logs/pipeline_locks' if lock_root is None else lock_root
 try:
  with (directory/f'{gid}_spincheck.lock').open('r') as stream:
   try:fcntl.flock(stream,fcntl.LOCK_SH|fcntl.LOCK_NB)
   except BlockingIOError:return True
 except FileNotFoundError:return False
 return False


def qualify_DFT_spin_gap(gap,matching_native):
 """Conservative publication gate, not a universal bound on thermal shifts."""
 if matching_native and abs(gap)<=.1:
  return 'spin_ensemble_required',False
 if not matching_native and abs(gap)<=.2:
  return 'matching_native_spin_thermals_required',False
 if gap < -1e-4:return 'wrong_ground_state',False
 return ('DFT_native_spin_checked' if matching_native else 'DFT_electronic_checked_thermal_uncertainty'),True


def helper_spin_report_admission(report,helper,chosen_identity=None):
 """Unfinished/stale helper calculations cannot establish the chosen spin."""
 if not report or report.get('status') not in ('spin_gap_verified','lower_spin_promoted'):
  return 'source_addressed_spin_report_pending',False
 if chosen_identity is not None and report.get('chosen_identity')!=chosen_identity:
  return 'source_addressed_spin_report_stale',False
 if helper.get('spin_validation_status') not in (None,'validated'):
  return 'alternate_spin_validation_pending',False
 if report.get('ensemble_primary_integration_pending'):
  return 'spin_ensemble_integration_pending',False
 return 'completed source-addressed spin comparison',True


def native_partition_dependencies_ok(report,chosen,alternative):
 """A thermal spin gap must bind both corrected partitions, not raw Hessians alone."""
 if not report.get('thermochemical_gap_validated'):return True
 import hashlib
 from pathlib import Path
 protocol=chosen.get('thermal_partition_protocol')
 if (not protocol or report.get('native_partition_protocol')!=protocol
     or alternative.get('thermal_partition_protocol')!=protocol):return False
 dependencies=report.get('native_partition_dependencies') or {}
 try:
  for role,record in [('chosen',chosen),('alternative',alternative)]:
   path=record.get('thermal_postprocessing_path')
   if not path or hashlib.sha256(Path(path).read_bytes()).hexdigest()!=dependencies.get(role):return False
 except OSError:return False
 return True




def state_spin_evidence(gid,state):
 """Completed source-bound spin evidence for one selected state; no calculations."""
 import hashlib,json,math
 from pathlib import Path
 from redox.core.common import UMA,DFT,read_result,apply_thermal_policy,free_energy
 from redox.core.protocol import load_record
 def safe_json(path):
  try:
   value=json.loads(Path(path).read_text())
   return value if isinstance(value,dict) else None
  except (OSError,ValueError):return None
 def finite(value):return type(value) in (float,int) and math.isfinite(value)
 def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
 uma=safe_json(UMA/gid/state/'result.json') or {}
 chosen=read_result(gid,state)
 from redox.validation.spin_ensemble import admitted_ensemble
 mixed=admitted_ensemble(chosen,DFT/gid/state)
 if mixed is not None:
  return dict(spin_confidence='DFT_native_equilibrium_ensemble',spin_gap_eV=mixed['gap_eV'],
              spin_gap_source='two qualified native spin surfaces; equilibrium partition',
              spin_evidence_path=chosen['spin_ensemble_policy_path'],spin_primary_admissible=True,
              single_surface_kinetics_validated=False)
 alternatives=[];pending=None
 if chosen and finite(chosen.get('e_smd_eV')):
  from redox.qm.dft import _opt_protocol_mismatch
  from redox.validation.spins import read_comparison,snapshot
  from redox.qm.native_thermal import cache_path
  report=read_comparison(gid,state)
  paths=list((DFT.parent/'spincheck'/gid).glob(state+'_m*/result.json'))
  for path in paths:
   alt=safe_json(path)
   if not alt or alt.get('charge')!=chosen.get('charge') or alt.get('mult')==chosen.get('mult'):continue
   reason,ok=helper_spin_report_admission(report,alt,snapshot(gid,state))
   if not ok:pending=(reason,str(path));continue
   if alt['mult']!=report.get('alternative_mult'):pending=('alternate_multiplicity_report_mismatch',str(path));continue
   if _opt_protocol_mismatch(alt,alt['charge']):pending=('alternate_optimizer_protocol_mismatch',str(path));continue
   sp=load_record(path.parent,alt['charge'],alt['mult'])
   if sp is None:pending=('alternate_active_SP_pending',str(path));continue
   if finite(sp.get('smd_spin_contam')) and abs(sp['smd_spin_contam'])>.1:
    pending=('alternate_spin_contamination',str(path));continue
   egap=sp['e_smd_eV']-chosen['e_smd_eV']
   if not finite(report.get('electronic_gap_eV')) or abs(egap-report['electronic_gap_eV'])>1e-8:
    pending=('alternate_energy_report_stale',str(path));continue
   admitted=apply_thermal_policy(dict(alt,e_smd_eV=sp['e_smd_eV']),path.parent)
   ga,qa=free_energy(admitted);gc,qc=free_energy(chosen)
   same=(admitted.get('thermal_definition')=='native_solution_dft_qrrho_spin'==chosen.get('thermal_definition')
          and bool(admitted.get('thermal_protocol_hash')) and bool(chosen.get('thermal_protocol_hash')))
   deps=report.get('native_dependencies') or {}
   if report.get('thermochemical_gap_validated'):
    if not native_partition_dependencies_ok(report,chosen,admitted):
     pending=('native_spin_partition_dependencies_stale',str(path));continue
    cp=cache_path(DFT/gid/state,chosen['charge'],chosen['mult'])[0]
    ap=cache_path(path.parent,alt['charge'],alt['mult'])[0]
    if not (cp.exists() and ap.exists() and deps.get('chosen')==sha(cp) and deps.get('alternative')==sha(ap)):
     pending=('native_spin_dependencies_stale',str(path));continue
   if same and qa==qc=='ok' and ga is not None and gc is not None:
    gap=ga-gc
    if not report.get('thermochemical_gap_validated') or not finite(report.get('thermochemical_gap_eV')) or abs(gap-report['thermochemical_gap_eV'])>1e-8:
     pending=('native_spin_report_pending_or_stale',str(path));continue
    alternatives.append((gap,'DFT matching native free energy',str(path)))
   else:alternatives.append((egap,'DFT active SMD electronic; thermal gap unresolved',str(path)))
 if alternatives:
  gap,method,path=min(alternatives,key=lambda x:x[0])
  confidence,admissible=qualify_DFT_spin_gap(gap,method=='DFT matching native free energy')
  return dict(spin_confidence=confidence,spin_gap_eV=gap,spin_gap_source=method,
              spin_evidence_path=path,spin_primary_admissible=admissible)
 if pending:
  return dict(spin_confidence=pending[0],spin_gap_eV=None,spin_gap_source='current completed spin report required',
              spin_evidence_path=pending[1],spin_primary_admissible=False)
 gap=uma.get('spin_gap_eV')
 if finite(gap):
  return dict(spin_confidence='DFT_check_required' if gap<.5 else 'UMA_large_gap_preliminary',
              spin_gap_eV=gap,spin_gap_source='UMA gas preoptimization; solution ground state not independently established',
              spin_evidence_path=str(UMA/gid/state/'result.json'),spin_primary_admissible=gap>=.5)
 standard=uma.get('mults_scanned')==[2] and uma.get('chosen_mult')==2
 return dict(spin_confidence='standard_doublet_model' if standard else 'spin_evidence_pending',
             spin_gap_eV=None,spin_gap_source='electron parity / no alternative scan' if standard else 'pending',
             spin_evidence_path=str(UMA/gid/state/'result.json') if uma else None,spin_primary_admissible=standard)

