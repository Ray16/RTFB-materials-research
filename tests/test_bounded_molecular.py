"""The shared launcher must bind a live budget owner before observation starts."""
import json
import os
from pathlib import Path

from redox.workflows import bounded_molecular, molecular_native_namespace, screening_control, response_diagnostic


def test_budget_observer_receives_owner_limits_and_frozen_sources(tmp_path, monkeypatch):
    case = tmp_path/'case'; case.mkdir()
    (case/'plan.json').write_text(json.dumps(dict(maximum_wall_seconds=14400, maximum_durable_output_bytes=2*2**30)))
    namespace, caller, out = tmp_path/'namespace.json', tmp_path/'caller.py', tmp_path/'placement'
    caller.write_text('# explicit caller source\n')
    spec = dict(case=str(case))
    events = []
    monkeypatch.setattr(molecular_native_namespace, 'verify', lambda p: spec)
    monkeypatch.setattr(screening_control, 'probe', lambda p: events.append('probe'))
    child = object()
    def execute(*args):
        events.append('execute'); out.mkdir()
        return child, dict(started=1)
    def observe(actual_child, actual_spec, result, output):
        events.append('observe')
        assert actual_child is child and actual_spec == spec and output == out
        assert result['budget_observer_owner']['pid'] == os.getpid()
        assert result['budget_observer_owner']['process_start'] == Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()[19]
        assert result['maximum_wall_seconds'] == 14400
        assert result['maximum_durable_output_bytes'] == 2*2**30
        assert json.loads((out/'launch.json').read_text()) == result
        assert len(result['bounded_launcher_sources']) == 3
        assert (out/'scientific_launcher_source.py').read_bytes() == caller.read_bytes()
        return dict(status='budget_handoff_verified')
    monkeypatch.setattr(molecular_native_namespace, 'launch', execute)
    monkeypatch.setattr(response_diagnostic, 'observe_budget', observe)
    assert bounded_molecular.launch(namespace, 0, out, caller_source=caller)['status'] == 'budget_handoff_verified'
    assert events == ['probe','execute','observe']
