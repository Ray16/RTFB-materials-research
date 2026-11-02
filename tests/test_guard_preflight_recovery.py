"""A timed-out observation cannot authorize recovery of a live/unknown owner."""
import json
from types import SimpleNamespace
import pytest
from redox.workflows import sampling_energy_recovery as recovery


def source(tmp_path,monkeypatch,live=False):
    root=tmp_path/'source';case=root/'case';guard=case/'residency_guard';placement=root/'placement_001'
    guard.mkdir(parents=True);placement.mkdir();namespace=root/'namespace.json';namespace.write_text('{}')
    owner=dict(host='test-host',pid=100,process_start='exact-start')
    (guard/'owner.json').write_text(json.dumps(owner));(case/'plan.json').write_text(json.dumps(dict(raw_root=str(case/'raw'))))
    (placement/'completion.json').write_text(json.dumps(dict(returncode=1,budget_stop=None)))
    (placement/'launch.json').write_text(json.dumps(dict(host='test-host',pid=101,process_start='wrapper-start',budget_observer_owner=dict(host='test-host',pid=102,process_start='controller-start'))))
    (placement/'command.log').write_text('TimeoutExpired: nvidia-smi initial guard observation')
    monkeypatch.setattr(recovery,'verify',lambda path:dict(case=str(case)))
    monkeypatch.setattr(recovery,'load',lambda *args:SimpleNamespace(alive=lambda owner:live))
    return namespace,case


def test_closed_guard_and_absent_science_are_explicitly_recorded(tmp_path,monkeypatch):
    namespace,_=source(tmp_path,monkeypatch)
    result=recovery.guard_preflight_source(namespace)
    assert result['checkpoint'] is None and not result['source_qualification_transferred']


@pytest.mark.parametrize('live',[True,None])
def test_live_or_unverified_process_never_restarts_after_observation_timeout(tmp_path,monkeypatch,live):
    namespace,_=source(tmp_path,monkeypatch,live)
    with pytest.raises(ValueError):recovery.guard_preflight_source(namespace)


@pytest.mark.parametrize('started',['status.json','host_guard','tensor_runtime','raw'])
def test_any_scientific_or_inner_runtime_start_requires_a_different_review(tmp_path,monkeypatch,started):
    namespace,case=source(tmp_path,monkeypatch)
    if started.endswith('.json'):(case/started).write_text('{}')
    else:(case/started).mkdir()
    with pytest.raises(ValueError):recovery.guard_preflight_source(namespace)
