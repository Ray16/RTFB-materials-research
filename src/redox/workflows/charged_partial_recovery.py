"""Recover a closed RKS H1 allocation failure from its exact partial checkpoint.

Only the saved solute partial Hessian and its original orbital/grid frame are
reused. Full solution H1, response, total Hessian and all physical checks remain
mandatory. The 16-GiB RKS workspace is a separate prospective H100 profile.
"""
import argparse
import json
import os
from pathlib import Path
import shutil

from redox.core.common import ROOT
from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json
from redox.workflows.molecular_native_namespace import verify,load
from redox.workflows.charged_polymer_pilot import prepare_native_component
from redox.workflows.native_component_resources import LARGE_CLASS,EXPANDED_RHF_CLASS
from redox.qm.hessian_restart import read_frame


def prepare(namespace,small_review,output):
    spec=verify(namespace);case=Path(spec['case']);old=json.loads((case/'plan.json').read_text())
    stage=json.loads((case/'status.json').read_text());guard=json.loads((case/'residency_guard/owner.json').read_text())
    placements=list(case.parent.glob('placement*/completion.json'))
    if len(placements)!=1:raise ValueError('One closed larger pilot attempt required')
    completion=json.loads(placements[0].read_text())
    owners=load(ROOT/'scripts/pipeline/observe_archived_native_namespace.py','partial_recovery_closed_owner')
    if (old.get('native_resource_class')!=LARGE_CLASS or old['charge']!=4 or old['mult']!=1
            or stage['status']!='charged_polymer_native_component_failed_preserved'
            or not stage.get('error','').startswith("MemoryError('Insufficient GPU memory.")
            or completion['returncode'] in (None,0) or completion['budget_stop'] is not None
            or not guard.get('guard_stopped') or owners.alive(stage['owner']) is not False
            or any(guard.get(k)!=stage['owner'].get(k) for k in ('host','pid','process_start'))
            or (case/'Hessian.npy').exists() or (case/'execution/rks_make_h1_state_identity.json').exists()):
        raise ValueError('Proven closed allocation failure before completed H1, with partial output preserved, required')
    frame,arrays=read_frame(case,partial_only=True)
    from redox.validation.scf_frame import checkpoint_frame,checkpoint_spin
    checkpoint=case/'SCF_evaluations/0000/0000/SCF.h5'
    actual=checkpoint_frame(checkpoint,charge=4,mult=1,basis='def2-svp',
        coordinates_Bohr=frame['identity']['coordinates_Bohr'],symbols=old['symbols'],maximum_error_Bohr=1e-12)
    spin=checkpoint_spin(checkpoint)
    result=prepare_native_component(Path(old['source_stationary_review']),output,small_review,EXPANDED_RHF_CLASS)
    current=json.loads((output/'namespace.json').read_text());target=Path(current['case']);sources=Path(current['sources'])
    plan=json.loads((target/'plan.json').read_text())
    if (plan['solvent_model']!=old['solvent_model'] or plan['numerical_settings']!=old['numerical_settings']
            or sha256_file(sources/'input.xyz')!=sha256_file(Path(spec['sources'])/'input.xyz')):
        raise ValueError('Partial recovery must retain the exact geometry, physical model and numerical settings')
    evidence=[namespace,case/'plan.json',case/'status.json',case/'residency_guard/owner.json',placements[0]]
    proof=dict(status='exact_partial_checkpoint_frame_verified_whole_Hessian_unqualified',
        source_namespace=str(namespace),source_sha256={str(p):sha256_file(p) for p in evidence},
        frame=frame,actual_checkpoint=actual,wavefunction_spin=spin,
        partial_only=True,full_solution_H1_recomputation_required=True,physical_gates_changed=False,
        native_G_qualified=False,material_admission=False)
    atomic_write_json(sources/'partial_checkpoint_review.json',proof)
    shutil.copyfile(Path(__file__),sources/'partial_recovery_preparer.py')
    plan.update(resume_partial_Hessian=dict(source_namespace=str(namespace),source_case=str(case),frame=frame),
        restart_scope='saved solute partial tensor only; original grids/orbitals; full solution H1/response recomputed',
        prior_allocation_failure_preserved=True)
    plan['source_sha256'].update(frame['source_sha256']);plan['source_sha256'].update(proof['source_sha256'])
    plan['source_sha256'].update({str(sources/name):sha256_file(sources/name) for name in ('partial_checkpoint_review.json','partial_recovery_preparer.py')})
    atomic_write_json(target/'plan.json',plan);current['plan_sha256']=sha256_file(target/'plan.json');atomic_write_json(output/'namespace.json',current)
    verify(output/'namespace.json')
    return dict(status='larger_native_partial_checkpoint_recovery_prepared',namespace=str(output/'namespace.json'),
        original=str(namespace),source_frame_verified=True,resource_profile=EXPANDED_RHF_CLASS,
        full_solution_H1_recomputation_required=True,physical_gates_changed=False,quantum_jobs_submitted=0)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('spec','small-pilot-review','out'):p.add_argument('--'+key,type=Path,required=True)
    a=p.parse_args()
    from redox.workflows.resources import allocation_guard,host_capacity,cpu_allocation
    with allocation_guard():host_capacity(1,minimum_ram_gib=4);os.sched_setaffinity(0,cpu_allocation(1))
    print(json.dumps(prepare(a.spec.resolve(),a.small_pilot_review.resolve(),a.out.resolve())))
