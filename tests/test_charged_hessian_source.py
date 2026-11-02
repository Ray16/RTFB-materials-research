"""A later I/O error cannot admit an incomplete response or bypass ownership/QC."""
from copy import deepcopy
from pathlib import Path

import pytest

from redox.validation.charged_hessian_source import response_calls,io_closure


def solver():
    return dict(status='response_solver_failed_preserved',equations_changed=False,level_shift_added=False,
        calls=[dict(shape=[3,9],solver_info=0,precision='float64',level_shift_added=False,
                    true_residual_norm=1e-8,solver_absolute_target=1e-7)]*3)


def test_complete_solved_rows_survive_only_an_explicit_later_failure_scope():
    response_calls(solver(),3,late_io_failure=True)
    with pytest.raises(ValueError):response_calls(solver(),3)


@pytest.mark.parametrize('change',['missing_row','info','nan_residual','infinite_target','loose_residual','shift','precision'])
def test_failed_or_missing_response_is_never_recovered_as_complete(change):
    data=deepcopy(solver());row=data['calls'][0]
    if change=='missing_row':data['calls'].pop()
    if change=='info':row['solver_info']=1
    if change=='nan_residual':row['true_residual_norm']=float('nan')
    if change=='infinite_target':row['solver_absolute_target']=float('inf')
    if change=='loose_residual':row['true_residual_norm']=1e-3
    if change=='shift':row['level_shift_added']=True
    if change=='precision':row['precision']='float32'
    with pytest.raises(ValueError):response_calls(data,3,late_io_failure=True)


@pytest.fixture
def closure():
    owner=dict(host='node',pid=11,process_start='123')
    error='BlockingIOError(11, "unable to lock file")'
    stage=dict(status='charged_polymer_native_component_failed_preserved',owner=owner,error=error,
        native_G_qualified=False,observed_peak_GPU_MiB=30000)
    launch=dict(host='node',CPU_affinity=[2,3])
    completion=dict(returncode=1,budget_stop=None)
    guard=dict(owner,CPU_affinity=[2,3],guard_stopped=True)
    residency=dict(own_residency_observed=True,foreign_resident_PIDs=[])
    host=dict(owner,CPU_affinity=[2,3],error=error,status='effective_memory_guard_failed_preserved',
        maximum_RSS_GiB=80,minimum_free_GiB=32,peak_tree_RSS_bytes=40*2**30,
        self_VmHWM_bytes=40*2**30,free_RAM_bytes=100*2**30)
    runtime=dict(pid=11,CPU_affinity=[2,3],error=error,status='runtime_or_task_failed_preserved',program='/worker.py')
    return [stage,launch,completion,guard,residency,host,runtime]


def test_exact_closed_io_failure_preserves_matrix_for_validation(closure):
    io_closure(*closure,owners_closed=True,worker=Path('/worker.py'),gpu_limit=88000)


@pytest.mark.parametrize('change',['live','foreign','budget','memory','wrong_error','wrong_start','wrong_worker','unclosed'])
def test_uncertain_ownership_resources_or_other_failures_remain_unqualified(closure,change):
    stage,launch,completion,guard,residency,host,runtime=closure
    closed=True
    if change=='live':closed=False
    if change=='foreign':residency['foreign_resident_PIDs']=[22]
    if change=='budget':completion['budget_stop']='deadline'
    if change=='memory':host['self_VmHWM_bytes']=90*2**30
    if change=='wrong_error':stage['error']='ValueError("failed minimum")'
    if change=='wrong_start':guard['process_start']='456'
    if change=='wrong_worker':runtime['program']='/other.py'
    if change=='unclosed':completion['returncode']=None
    with pytest.raises(ValueError):io_closure(*closure,owners_closed=closed,worker=Path('/worker.py'),gpu_limit=88000)
