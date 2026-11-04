"""Only a closed physical failure can authorize corrective geometry work."""
from pathlib import Path
import pytest
from redox.validation.failed_execution import physical_failure


def records(tmp_path):
    case=tmp_path/'case';case.mkdir();placement=tmp_path/'placement_001';placement.mkdir();(placement/'launch.json').write_text('{}')
    identity=dict(pid=11,host='node',process_start='start');cpus=[2,3]
    data={case/'status.json':dict(status='failed',error='physical reason',native_G_qualified=False,owner=identity),
        placement/'launch.json':dict(host='node',pid=12,process_start='wrapper',budget_observer_owner=dict(host='node',pid=13,process_start='controller'),CPU_affinity=cpus),
        placement/'completion.json':dict(returncode=1,budget_stop=None),
        case/'residency_guard/owner.json':dict(identity,guard_stopped=True,CPU_affinity=cpus),
        case/'residency_guard/state.json':dict(foreign_resident_PIDs=[],own_residency_observed=True),
        case/'host_guard/status.json':dict(identity,CPU_affinity=cpus,maximum_RSS_GiB=40,minimum_free_GiB=32,peak_tree_RSS_bytes=2**30,self_VmHWM_bytes=2**30,free_RAM_bytes=64*2**30),
        case/'tensor_runtime/status.json':dict(pid=11,status='runtime_or_task_failed_preserved',program=str(tmp_path/'sources/worker.py'),CPU_affinity=cpus)}
    return case,data


def test_closed_actual_failure_is_retained_as_evidence(tmp_path):
    case,data=records(tmp_path)
    assert physical_failure(case,data.__getitem__,lambda owner:False,status='failed',error='physical reason',host_GiB=40)['owner']['pid']==11


@pytest.mark.parametrize('bad',['live','unknown','budget','foreign','runtime','owner','memory','error'])
def test_resource_owner_or_scope_uncertainty_cannot_be_promoted_to_a_physical_failure(tmp_path,bad):
    case,data=records(tmp_path);alive=False
    if bad=='live':alive=True
    if bad=='unknown':alive=None
    if bad=='budget':data[tmp_path/'placement_001/completion.json']['budget_stop']='wall_time_budget'
    if bad=='foreign':data[case/'residency_guard/state.json']['foreign_resident_PIDs']=[99]
    if bad=='runtime':data[case/'tensor_runtime/status.json']['status']='running'
    if bad=='owner':data[case/'host_guard/status.json']['pid']=100
    if bad=='memory':data[case/'host_guard/status.json']['peak_tree_RSS_bytes']=41*2**30
    if bad=='error':data[case/'status.json']['error']='different failure'
    with pytest.raises(ValueError):physical_failure(case,data.__getitem__,lambda owner:alive,status='failed',error='physical reason',host_GiB=40)
