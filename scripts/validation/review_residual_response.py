#!/usr/bin/env python
"""Bind full-spin Hessian replays and a fixed-state large response diagnostic."""
import argparse
import json
from pathlib import Path
import time

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api


def checked_controls(root):
    path=Path(root)/'comparison.json';row=json.loads(path.read_text())
    if (row['status']!='full_spin_replays_and_large_first_atom_response_checked'
            or row['full_native_G_qualified'] or row['equations_changed']):
        raise ValueError('Full spin and actual large first-atom checks are required')
    api.verify_sources(row['source_sha256'])
    return {str(path):api.digest(path),**row['source_sha256']}


def review(root,probe,baseline):
    import numpy as np
    small=api.load_module('independent_residual_spin_checks',ROOT/'scripts/validation/probe_fixed_orbital_response.py')
    proof=small.controls(root)
    original=api.load_module('retained_bounded_DF_and_cross_engine_checks',ROOT/'scripts/validation/review_cutensor_controls.py')
    proof.update(original.checked_controls(baseline))
    path=probe/'response/status.json';row=json.loads(path.read_text())
    if (row['status']!='first_atom_response_equations_checked_full_native_QC_pending'
            or not row.get('finished') or row['SCF_recomputed'] or row['full_Hessian_qualified'] or row['native_G_qualified']):
        raise ValueError('Completed isolated response diagnostic required; never a full Hessian/G proxy')
    owner=api.load_module('response_diagnostic_ended_owner',ROOT/'scripts/pipeline/run_polymer_spacing.py')
    if owner.owner_alive(row['owner']):raise ValueError('Probe owner still active')
    api.verify_sources(row['output_sha256'])
    plan=api.read_plan(probe);api.verify_sources(plan['producer_sha256']);api.verify_sources(plan['input_sha256'])
    raw=Path(plan['raw_root']);residual=json.loads((raw/'physical_residual.json').read_text())
    solves=json.loads((raw/'execution/response_solver.json').read_text())
    if (not plan.get('equilibrium_solvation') or not plan.get('PCM_nuclear_H1_required')
            or not residual['passed'] or not np.isfinite(residual['maximum_Ha']) or residual['maximum_Ha']>1e-8
            or residual['all_atoms_solved'] or not solves['calls']
            or any(c['solver_info']!=0 or c['true_residual_norm']>c['solver_absolute_target'] for c in solves['calls'])
            or not 0<row['observed_peak_GPU_MiB']<=31744):
        raise ValueError('First-atom physical equations and measured resources must pass')
    runtime_path=probe/'tensor_runtime/status.json';runtime=json.loads(runtime_path.read_text())
    if runtime['status']!='task_returned_under_checked_cutensor_runtime' or runtime['pid']!=row['owner']['pid']:
        raise ValueError('Actual completed contraction runtime differs from the probe')
    guard=probe/'residency_guard';r=json.loads((guard/'state.json').read_text())
    if r['foreign_resident_PIDs'] or not r['own_residency_observed'] or (guard/'interruption.json').exists():
        raise ValueError('Large response probe did not retain exclusive residency')
    proof.update(row['output_sha256']);proof.update(plan['producer_sha256']);proof.update(plan['input_sha256'])
    proof.update({str(p):api.digest(p) for p in (Path(__file__).resolve(),path,probe/'plan.json',runtime_path,
        guard/'owner.json',guard/'state.json')})
    controls=[]
    for m in (1,2):
        row=json.loads((root/'controls'/str(m)/'status.json').read_text())
        controls.append(dict(mult=m,differences=row['differences'],peak_GPU_MiB=row['observed_peak_GPU_MiB']))
    result=dict(status='full_spin_replays_and_large_first_atom_response_checked',checked=time.time(),
        controls=controls,large_first_atom_physical_residual_Ha=residual['maximum_Ha'],
        runtime_profile=plan['runtime_profile'],baseline_controls=str(baseline),source_sha256=proof,
        equations_changed=False,level_shift_added=False,full_native_G_qualified=False,
        scope='Full two-spin implementation equivalence and one exact large response; full molecular native QC remains required')
    atomic_write_json(root/'comparison.json',result)
    print(json.dumps(dict(status=result['status'],first_atom_residual_Ha=residual['maximum_Ha'])))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ('root','probe','baseline'):parser.add_argument('--'+name,type=Path,required=True)
    args=parser.parse_args()
    review(args.root.resolve(),args.probe.resolve(),args.baseline.resolve())


if __name__=='__main__':
    main()
