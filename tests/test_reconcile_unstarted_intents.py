import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace

import pytest


source = Path(__file__).resolve().parents[1] / 'scripts/pipeline/reconcile_unstarted_intents.py'
spec = importlib.util.spec_from_file_location('reconcile_intents', source)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


@pytest.mark.parametrize('marker', ['launch.json', 'dispatch_attempts', 'experimental_energy',
                                    'residency_guard', 'host_guard', 'tensor_runtime', 'raw'])
def test_any_attempt_artifact_prevents_reconciliation(tmp_path, marker):
    case = tmp_path / 'case'
    case.mkdir()
    raw = tmp_path / 'raw'
    assert module.no_attempt(case, raw)
    (raw if marker == 'raw' else case / marker).mkdir()
    assert not module.no_attempt(case, raw)


def test_remote_query_transmits_case_paths_only_via_stdin(monkeypatch):
    def fake_run(argv, **kwargs):
        assert '/sensitive/case' not in str(argv)
        payload = json.loads(kwargs['input'])
        assert payload['cases']['test'] == '/sensitive/case'
        assert kwargs['check'] and kwargs['timeout'] == 15
        return SimpleNamespace(stdout='{"matches": [], "old_pid_starts": {}}')
    monkeypatch.setattr(module.subprocess, 'run', fake_run)
    assert module.remote_check('lambda1', {'test': '/sensitive/case'}, [12])['matches'] == []
