from copy import deepcopy
from pathlib import Path

import pytest

from redox.validation.execution_review import validate_closed


@pytest.fixture
def records():
    owner=dict(host='test-host',pid=10,process_start='100')
    cpus=[2,3]
    return [dict(status='test_finished',native_G_qualified=False,observed_peak_GPU_MiB=23000,owner=owner),
        dict(host='test-host',CPU_affinity=cpus),dict(returncode=0,budget_stop=None),
        dict(**owner,guard_stopped=True,CPU_affinity=cpus),
        dict(own_residency_observed=True,foreign_resident_PIDs=[]),
        dict(**owner,status='effective_memory_guard_finished_profile_review_pending',CPU_affinity=cpus,
            maximum_RSS_GiB=40,minimum_free_GiB=32,peak_tree_RSS_bytes=2**30,self_VmHWM_bytes=2**30,free_RAM_bytes=100*2**30),
        dict(pid=10,status='task_returned_under_checked_cutensor_runtime',CPU_affinity=cpus,program='/test/worker.py')]


def validate(rows,alive=False):
    validate_closed(*rows,owner_alive=alive,expected_status='test_finished',host_GiB=40,GPU_MiB=31744,worker=Path('/test/worker.py'))


def test_normal_exact_execution_closure_is_accepted(records):validate(records)


@pytest.mark.parametrize('index,field,value',[
    (0,'observed_peak_GPU_MiB',float('nan')),(0,'observed_peak_GPU_MiB',40000),
    (0,'native_G_qualified',True),(2,'returncode',1),(2,'budget_stop','wall_time_budget'),
    (3,'process_start','101'),(3,'host','other-host'),(3,'guard_stopped',False),
    (4,'foreign_resident_PIDs',[20]),(5,'free_RAM_bytes',float('nan')),
    (5,'peak_tree_RSS_bytes',41*2**30),(5,'CPU_affinity',[2,2]),
    (6,'program','/test/other.py'),(6,'pid',20),
])
def test_broken_identity_nonfinite_resources_or_abnormal_exit_cannot_pass(records,index,field,value):
    rows=deepcopy(records);rows[index][field]=value
    with pytest.raises(ValueError):validate(rows)


def test_live_or_unknown_owner_cannot_be_accepted(records):
    for alive in (True,None):
        with pytest.raises(ValueError):validate(records,alive)


def test_source_capture_keeps_as_used_bytes_after_the_source_changes(tmp_path):
    from redox.validation.execution_review import freeze_source
    source=tmp_path/'review.py';source.write_bytes(b'first implementation')
    captured=source.read_bytes();source.write_bytes(b'new implementation')
    first=freeze_source(tmp_path/'case',source,captured)
    current=freeze_source(tmp_path/'case',source)
    assert first.read_bytes()==captured and current.read_bytes()==b'new implementation'
    assert first!=current
