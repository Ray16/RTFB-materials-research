import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace
import pytest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('geometry_intents', ROOT / 'scripts/pipeline/reconcile_geometry_intents.py')
mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)


def test_scientific_attempts_never_reopened(tmp_path):
    case = tmp_path / 'case'; case.mkdir(); raw = tmp_path / 'raw'
    assert mod.no_attempt(case, raw)
    for name in mod.FORBIDDEN:
        p = case / name; p.mkdir()
        assert not mod.no_attempt(case, raw)
        p.rmdir()
    raw.mkdir()
    assert not mod.no_attempt(case, raw)


def test_remote_unknown_or_matching_process_fails_closed(monkeypatch):
    for row in [dict(host='x', matches=[], errors=[{'pid':1}]),
                dict(host='x', matches=[{'pid':2}], errors=[]),
                dict(host='wrong', matches=[], errors=[])]:
        monkeypatch.setattr(mod.subprocess, 'run', lambda *a, row=row, **k: SimpleNamespace(stdout=json.dumps(row)))
        with pytest.raises(ValueError):mod.remote_check('x', {'a':['/case']}, [1])


@pytest.mark.parametrize('block', ['none', 'live_owner', 'live_gate'])
def test_reconciliation_preserves_original_before_reopening(tmp_path, monkeypatch, block):
    case = tmp_path / 'case'; case.mkdir(); raw = tmp_path / 'raw'
    original = b'{"host":"x","GPU":0,"owner":{"host":"x","pid":10,"process_start":"old"}}\n'
    (case / 'dispatch_intent.json').write_bytes(original); (case / 'plan.json').write_text('{}')
    audit = tmp_path / 'input.json'; audit.write_text('{}'); gate = tmp_path / 'gate'; gate.write_text('gate'); policy = tmp_path / 'policy'; policy.write_text('policy')
    inventory = SimpleNamespace(GATE=gate, POLICY=policy, permitted_hosts=lambda p:{'x'})
    monkeypatch.setattr(mod.api, 'load_module', lambda *a:inventory)
    monkeypatch.setattr(mod, 'collect', lambda a:{'a':dict(case=case, raw=raw, intent_bytes=original, intent=json.loads(original), plan_sha256=mod.api.digest(case/'plan.json'))})
    monkeypatch.setattr(mod, 'remote_check', lambda *a:dict(host='x', matches=[], errors=[], old_pid_starts={'10':'old'} if block == 'live_owner' else {}))
    monkeypatch.setattr(mod.subprocess, 'run', lambda *a, **k:SimpleNamespace(stdout='x_gpu0: 10 x user label [live]' if block == 'live_gate' else ''))
    out = tmp_path / 'review'
    if block != 'none':
        with pytest.raises(ValueError):mod.reconcile(audit, out)
        assert (case/'dispatch_intent.json').read_bytes() == original and not out.exists()
        return
    result = mod.reconcile(audit, out)
    assert result['reopened'] == 1 and (out/'a.intent.json').read_bytes() == original
    assert json.loads((case/'dispatch_intent.json').read_text())['outcome'] == 'not_started'
    assert not raw.exists()
