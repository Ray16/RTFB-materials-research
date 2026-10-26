import importlib.util
from pathlib import Path

import pytest


source = Path(__file__).resolve().parents[1] / 'scripts/validation/wait_uma_pair_cpu_preflight.py'
spec = importlib.util.spec_from_file_location('uma_pair_cpu_preflight', source)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def test_only_canonical_helper_is_replaced():
    command = ['python', str(module.ORIGINAL), 'lambda13', '--module', 'torch']
    revised = module.preflight_command(command)
    assert revised == ['python', str(module.CPU_ONLY), 'lambda13', '--module', 'torch']
    assert command[1] == str(module.ORIGINAL)
    assert module.preflight_command(['ssh', 'lambda13', 'true']) == ['ssh', 'lambda13', 'true']
    with pytest.raises(ValueError, match='Ambiguous'):
        module.preflight_command([str(module.ORIGINAL), str(module.ORIGINAL)])


def test_process_local_patch_is_restored_on_failure(monkeypatch):
    original_verify, original_run = module.api.verify_sources, module.subprocess.run
    def fail(*args, **kwargs):
        assert module.api.verify_sources is not original_verify
        assert module.subprocess.run is not original_run
        raise RuntimeError('synthetic waiter failure')
    monkeypatch.setattr(module.runpy, 'run_path', fail)
    with pytest.raises(RuntimeError, match='synthetic'):
        module.main()
    assert module.api.verify_sources is original_verify
    assert module.subprocess.run is original_run
