import importlib.util
import json
from pathlib import Path
import subprocess
from types import SimpleNamespace

import pytest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('measured_energy_preparation_transport',
                                             ROOT/'scripts/pipeline/prepare_measured_energy_class.py')
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


@pytest.fixture(autouse=True)
def clear_stopped_owner_cache():
    module.STOPPED_OWNERS.clear()


def remote_owner():
    return dict(host='a-remote-node', pid=12345, process_start='fixed-start')


@pytest.mark.parametrize('reply,expected', [('True\n', True), ('False\n', False)])
def test_owner_query_preserves_exact_identity_with_portable_remote_python(monkeypatch, reply, expected):
    calls = []

    def run(argv, **kwargs):
        calls.append(argv)
        return SimpleNamespace(stdout=reply)

    monkeypatch.setattr(module.subprocess, 'run', run)
    assert module.owner_alive(remote_owner()) is expected
    assert '/usr/bin/python3 -I -S' in calls[0][-1]
    assert calls[0][2:4] == ['-S', 'none']
    assert '/proc/12345/stat' in calls[0][-1] and 'fixed-start' in calls[0][-1]


@pytest.mark.parametrize('reply', ['unreachable', ''])
def test_unverified_owner_cannot_be_treated_as_stopped(monkeypatch, reply):
    monkeypatch.setattr(module.subprocess, 'run', lambda *a, **k: SimpleNamespace(stdout=reply))
    with pytest.raises(module.OwnershipUnverified):
        module.owner_alive(remote_owner())


def test_rpc_timeout_preserves_unknown_ownership(monkeypatch):
    def fail(*args, **kwargs):
        raise subprocess.TimeoutExpired('ssh', 20)

    monkeypatch.setattr(module.subprocess, 'run', fail)
    with pytest.raises(module.OwnershipUnverified):
        module.owner_alive(remote_owner())


def test_transport_adapter_restores_original_loader_after_failure(monkeypatch):
    original = module.api.load_module
    with pytest.raises(RuntimeError):
        with module.ownership_transport():
            assert module.api.load_module is not original
            raise RuntimeError('preparation interrupted')
    assert module.api.load_module is original


def test_geometry_recovery_cannot_change_the_declared_task(tmp_path, monkeypatch):
    batch = tmp_path/'source'
    batch.mkdir()
    (batch/'batch.json').write_text(json.dumps(dict(entries=[dict(key='k', campaign='/recovered')],
                                                  source_sha256={'bound': 'source'})))
    monkeypatch.setattr(module.api, 'verify_sources', lambda bindings: None)
    monkeypatch.setattr(module.api, 'read_plan', lambda path: dict(tasks=[dict(charge=0 if str(path)=='/old' else 1)]))
    cfg = dict(geometry_subset_inputs=[dict(batch=str(batch), keys=['k'])], recovered_from=dict(k='/old'))
    with pytest.raises(ValueError, match='Recovery changed'):
        module.prepare_recovered(cfg, tmp_path/'configuration.json', None, tmp_path/'out')
