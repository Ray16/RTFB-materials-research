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


