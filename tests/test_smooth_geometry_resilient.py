import importlib.util
from pathlib import Path
from types import SimpleNamespace

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('smooth_geometry_resilient', ROOT / 'scripts/pipeline/coordinate_smooth_geometry_resilient.py')
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


def test_cache_requires_matching_policy_hosts_and_recent_time():
    row = dict(hosts=['a', 'b'], policy_sha256='x', checked=100, slots=[], checks=[])
    assert mod.cache_valid(row, ['b', 'a'], 'x', 130)
    assert not mod.cache_valid(row, ['a'], 'x', 110)
    assert not mod.cache_valid(row, ['a', 'b'], 'y', 110)
    assert not mod.cache_valid(row, ['a', 'b'], 'x', 131)
    assert not mod.cache_valid(row, ['a', 'b'], 'x', 99)


def test_inventory_preserves_unverified_hosts_and_reuses_cache(tmp_path, monkeypatch):
    monkeypatch.setattr(mod, 'CACHE', tmp_path / 'cache')
    policy = tmp_path / 'policy'; policy.write_text('policy')
    calls = []
    def query(hosts):
        calls.append(hosts)
        return [('a', 1)], [dict(host='a', status='verified'), dict(host='b', status='inventory_unverified')]
    inv = SimpleNamespace(POLICY=policy, query_all=query)
    slots, checks = mod.shared_inventory(['a', 'b'], inv)
    assert slots == [('a', 1)] and checks[1]['status'] == 'inventory_unverified'
    assert mod.shared_inventory(['b', 'a'], inv)[0] == slots
    assert len(calls) == 1
    policy.write_text('different')
    mod.shared_inventory(['a', 'b'], inv)
    assert len(calls) == 2


def test_adapter_rewrites_only_expected_existing_queue():
    module = mod.api.load_module('test_source_geometry_queue', mod.ORIGINAL)
    profile = {'library_path': '/example', 'python': '/example/python'}
    mod.adapt(module, profile)
    assert module.SCRIPT == Path(mod.__file__).resolve()
    assert callable(module.preflight) and callable(module.launch_local) and callable(module.coordinate)
