"""Evidence for retaining a completed ionic Hessian after later file-lock failure.

This admits a matrix for further validation, never a minimum or free energy.
The original failure, response frame and resource/ownership chain stay intact.
"""
import json
import math
from pathlib import Path

import numpy as np

from redox.core.hashing import sha256_file
from redox.validation.execution_review import host_usage_peak
from redox.workflows.molecular_native_namespace import verify
from redox.workflows.native_component_resources import peak_limit


def response_calls(solver, count, *, late_io_failure=False):
    expected=('response_solver_failed_preserved' if late_io_failure else
              'response_equations_solved_full_physical_QC_separate')
    calls=solver['calls']
    if (solver['status']!=expected or solver['equations_changed'] or solver['level_shift_added']
            or not calls or sum(row['shape'][0] for row in calls)!=3*count
            or any(row['solver_info']!=0 or row['precision']!='float64' or row['level_shift_added']
                   or type(row['shape'][0]) is not int or row['shape'][0]<=0
                   or not math.isfinite(row['true_residual_norm'])
                   or not math.isfinite(row['solver_absolute_target'])
                   or not 0<=row['true_residual_norm']<=row['solver_absolute_target']
                   or row['solver_absolute_target']<=0 for row in calls)):
        raise ValueError('All original-equation response rows must have solved before matrix reuse')


def full_response(case, plan, read, used, *, late_io_failure=False):
    case=Path(case);count=len(plan['symbols'])
    solver=read(case/'execution/response_solver.json')
    response_calls(solver,count,late_io_failure=late_io_failure)
    solution=read(case/'solution_H1/identity.json');path=Path(solution['path']);used.append(path)
    center=case/'points/0000/evaluation.json';grids=case/'grids.npz';used.extend([center,grids])
    if (solution['scope']!='full SMD H1 including PCM contribution' or not solution['equilibrium_solvation']
            or solution['file_sha256']!=sha256_file(path)
            or path.resolve()!=(case/'solution_H1/full_H1.npy').resolve()
            or solution['identity']['center_evaluation_sha256']!=sha256_file(center)
            or solution['identity']['grids_sha256']!=sha256_file(grids)
            or solution['identity']['screening_variant']!='corrected'
            or solution['identity']['source_namespace']!=plan['original_namespace']
            or solution['shape']!=[count,3,plan['response_nao'],plan['response_nocc']]):
        raise ValueError('Complete solvent H1 must retain its exact source orbital/grid frame')
    return dict(response_calls=len(solver['calls']),response_rows=3*count,
                source_response_independently_checked=True,full_response_recomputed=False)


def io_closure(stage, launch, completion, guard, residency, host, runtime, *,
               owners_closed, worker, gpu_limit):
    error=stage.get('error','');owner=stage['owner']
    if (not error.startswith('BlockingIOError(11,') or 'unable to lock file' not in error
            or host.get('error')!=error or runtime.get('error')!=error
            or stage['status']!='charged_polymer_native_component_failed_preserved'
            or stage['native_G_qualified'] is not False or owners_closed is not True
            or completion.get('returncode') in (None,0) or completion.get('budget_stop') is not None
            or not guard.get('guard_stopped') or not residency.get('own_residency_observed')
            or residency.get('foreign_resident_PIDs')
            or host.get('status')!='effective_memory_guard_failed_preserved'
            or runtime.get('status')!='runtime_or_task_failed_preserved'
            or runtime.get('program')!=str(worker)
            or host.get('maximum_RSS_GiB')!=80 or host.get('minimum_free_GiB')!=32
            or not math.isfinite(stage['observed_peak_GPU_MiB'])
            or not 0<stage['observed_peak_GPU_MiB']<=gpu_limit
            or host_usage_peak(host)>80*2**30
            or not math.isfinite(host['free_RAM_bytes']) or host['free_RAM_bytes']<32*2**30):
        raise ValueError('Only a closed, uncontended, within-budget checkpoint-file lock failure is recoverable')
    if (not owner['pid']==guard['pid']==host['pid']==runtime['pid']
            or not owner['process_start']==guard['process_start']==host['process_start']
            or not owner['host']==guard['host']==host['host']==launch['host']
            or not guard['CPU_affinity']==host['CPU_affinity']==runtime['CPU_affinity']==launch['CPU_affinity']
            or len(guard['CPU_affinity'])!=2 or len(set(guard['CPU_affinity']))!=2):
        raise ValueError('Exact source worker and two-CPU ownership chain required')


def source_evidence(namespace, alive):
    namespace=Path(namespace);spec=verify(namespace);case=Path(spec['case']);used=[namespace]
    def read(path):
        path=Path(path);used.append(path);return json.loads(path.read_text())
    plan=read(case/'plan.json');stage=read(case/'status.json')
    if (not plan.get('charged_polymer_native_component') or plan['charge']!=4 or plan['mult']!=1
            or plan.get('saved_full_Hessian') or plan.get('Hessian_validation_continuation')):
        raise ValueError('Original full double-graft ionic Hessian source required')
    placements=list(case.parent.glob('placement*/launch.json'))
    if len(placements)!=1 or (case/'residency_guard/external_interruption.json').exists():
        raise ValueError('One uncontended source placement required')
    launch=read(placements[0]);completion=read(placements[0].parent/'completion.json')
    guard=read(case/'residency_guard/owner.json');residency=read(case/'residency_guard/state.json')
    host=read(case/'host_guard/status.json');runtime=read(case/'tensor_runtime/status.json')
    owners=[stage['owner'],dict(host=launch['host'],pid=launch['pid'],process_start=launch['process_start']),
            completion['budget_observer_owner']]
    io_closure(stage,launch,completion,guard,residency,host,runtime,
        owners_closed=all(alive(owner) is False for owner in owners),
        worker=Path(spec['sources'])/'worker.py',gpu_limit=peak_limit(plan))
    matrix_path=case/'Hessian.npy';used.append(matrix_path)
    matrix=np.load(matrix_path,allow_pickle=False);count=len(plan['symbols'])
    center=read(case/'points/0000/evaluation.json');diagnostics=read(case/'Hessian_diagnostics.json')
    gradient=np.asarray(center['gradient_Ha_Bohr'])
    if (matrix.dtype!=np.float64 or matrix.shape!=(3*count,3*count) or not np.isfinite(matrix).all()
            or abs(matrix-matrix.T).max()>1e-6 or not diagnostics['all_finite']
            or diagnostics['maximum_asymmetry']>1e-6 or not center['SCF_converged']
            or gradient.shape!=(count,3) or not np.isfinite(gradient).all() or abs(gradient).max()>.00045):
        raise ValueError('Complete finite symmetric Hessian and stationary center required')
    proof=full_response(case,plan,read,used,late_io_failure=True)
    points=sorted((case/'points').glob('*/evaluation.json'));columns=read(case/'columns.json')
    if not 1<len(points)<1+4*len(columns['columns']):
        raise ValueError('Failure must leave incomplete post-Hessian derivative validation')
    used.extend(points)
    proof.update(status='closed_ionic_Hessian_retained_for_fresh_validation',source_namespace=str(namespace),
        source_case=str(case),source_charge=plan['charge'],source_mult=plan['mult'],
        saved_completed_evaluations=len(points),fresh_all_sampled_gradients_required=True,
        native_G_qualified=False,source_sha256={str(p):sha256_file(p) for p in used})
    verify(namespace)
    return spec,plan,proof


def saved_response(case, plan, read, used, alive):
    saved=plan['saved_full_Hessian']
    _,old,proof=source_evidence(Path(saved['namespace']),alive)
    used.extend(Path(p) for p in proof['source_sha256'])
    if (proof['source_sha256']!=saved['source_sha256'] or saved.get('charge')!=plan['charge']
            or saved.get('mult')!=plan['mult'] or not saved.get('source_response_independently_checked')
            or any(plan[k]!=old[k] for k in ('symbols','charge','mult','numerical_settings','solvent_model','prepared_case'))):
        raise ValueError('Saved ionic response differs from its proven source state')
    if sha256_file(Path(case)/'Hessian.npy')!=sha256_file(saved['matrix']):
        raise ValueError('Validation must retain the exact complete source matrix')
    reuse=read(Path(case)/'saved_matrix_reuse.json');calls=reuse['calls']
    if (reuse['source']!=saved or len(calls)!=1 or calls[0]['full_response_recomputed'] is not False
            or calls[0]['source_matrix_sha256']!=sha256_file(saved['matrix'])
            or not 0<=calls[0]['center_energy_difference_Ha']<=1e-8
            or not 0<=calls[0]['center_gradient_difference_Ha_Bohr']<=1e-6):
        raise ValueError('One exact matrix reuse and reproduced fresh center required')
    solver=read(Path(case)/'execution/response_solver.json')
    if (solver['status']!='no_new_response_solve_saved_full_Hessian_validation' or solver['calls']
            or solver['equations_changed'] or solver['level_shift_added'] or solver['source_response']!=saved['solver']):
        raise ValueError('Saved-matrix validation may not claim a new response solve')
    return proof
