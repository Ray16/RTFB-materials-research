import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace


SOURCE = Path(__file__).resolve().parents[1] / 'scripts/pipeline/coordinate_electronic_inventory_resilient.py'
spec = importlib.util.spec_from_file_location('resilient_electronic_inventory', SOURCE)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def test_unknown_host_stays_unknown_and_inventory_is_restored(tmp_path, monkeypatch):
    root = tmp_path / 'project'
    batch = root / 'batch'
    batch.mkdir(parents=True)
    (batch / 'batch.json').write_text(json.dumps({'source_sha256': {}}))
    profile = root / 'environments/cutensor22/runtime.json'
    profile.parent.mkdir(parents=True)
    profile.write_text(json.dumps({'python': module.sys.executable}))
    policy, gate = root / 'policy', root / 'gate'
    inventory = SimpleNamespace(__file__=str(root / 'inventory.py'), POLICY=policy, GATE=gate,
        permitted_hosts=lambda p: {'lambda1', 'lambda2'},
        query_all=lambda hosts: ([('lambda1', 2)], [
            {'host': 'lambda1', 'status': 'verified', 'slots': [('lambda1', 2)]},
            {'host': 'lambda2', 'status': 'inventory_unverified', 'slots': []}]))
    monkeypatch.setattr(module, 'ROOT', root)
    monkeypatch.setattr(module.api, 'load_module', lambda *a: inventory)
    monkeypatch.setattr(module.api, 'verify_sources', lambda x: None)
    monkeypatch.setattr(module.api, 'digest', lambda x: 'hash')
    monkeypatch.delenv('LD_LIBRARY_PATH', raising=False)
    monkeypatch.setattr(module.sys, 'argv', ['controller', '--batch', str(batch), '--hosts', 'lambda1,lambda2'])
    original = module.resources.free_slots
    def fake_controller(*args, **kwargs):
        assert module.resources.free_slots() == [('lambda1', 2)]
        status = json.loads((batch / 'resilient_inventory/status.json').read_text())
        assert status['unverified_hosts'] == ['lambda2']
        assert status['available_slots'] == [['lambda1', 2]]
    monkeypatch.setattr(module.runpy, 'run_path', fake_controller)
    module.main()
    assert module.resources.free_slots is original
    owner = next((batch / 'resilient_inventory').glob('owner_*.json'))
    assert json.loads(owner.read_text())['status'] == 'electronic_resilient_inventory_adapter_stopped'
