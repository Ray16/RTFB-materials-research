import copy
import json

import pytest

from redox.workflows.archived_native_namespace import mount_command, verify_spec


def spec():
    return dict(prior_root='/workspace/redox', archive='/archive/redox',
                original_campaign='/workspace/redox/results/polymer_models/campaigns/batch/key',
                original_raw_parent='/workspace/redox/calcs/polymer_models/batch',
                attempt='/durable/new/campaign', raw_parent='/durable/new/raw')


def test_archive_and_system_are_read_only_with_two_fresh_write_mounts():
    command = mount_command(spec())
    assert command[:5] == ['bwrap', '--die-with-parent', '--ro-bind', '/', '/']
    assert command.count('--bind') == 2
    assert command.count('--ro-bind') == 2
    assert '--unshare-pid' not in command
    assert '--tmpfs' in command


@pytest.mark.parametrize('field,value', [
    ('original_campaign', '/workspace/redox/src'),
    ('original_raw_parent', '/workspace/redox/calcs/dft'),
    ('attempt', '/archive/redox/changed'),
    ('raw_parent', '/workspace/redox/changed'),
])
def test_refuses_write_mounts_over_scientific_history_or_code(field, value):
    value_spec = spec()
    value_spec[field] = value
    with pytest.raises(ValueError):
        mount_command(value_spec)


def test_preserved_plan_changes_block_execution_before_gpu_access(tmp_path):
    from redox.core.hashing import sha256_file
    attempt = tmp_path/'attempt'
    attempt.mkdir()
    (attempt/'plan.json').write_text('{}')
    (attempt/'input.xyz').write_text('input')
    value = dict(attempt=str(attempt), plan_sha256=sha256_file(attempt/'plan.json'),
                 geometry_sha256=sha256_file(attempt/'input.xyz'))
    path = tmp_path/'namespace.json'
    path.write_text(json.dumps(value))
    assert verify_spec(path) == value
    (attempt/'plan.json').write_text('{"changed":true}')
    with pytest.raises(ValueError, match='plan'):
        verify_spec(path)


def test_launcher_stays_alive_until_namespace_child_finishes(tmp_path, monkeypatch):
    import sys
    from redox.workflows import archived_native_namespace as module
    observed = []
    monkeypatch.setattr(sys, 'argv', ['native-namespace', 'launch', '--spec', str(tmp_path/'spec.json'),
                                    '--gpu', '0', '--out', str(tmp_path)])
    monkeypatch.setattr(module, 'launch_reserved', lambda *args: dict(pid=1234))
    monkeypatch.setattr(module.os, 'waitpid', lambda pid, flags: observed.append((pid, flags)) or (pid, 0))
    module.main()
    assert observed == [(1234, 0)]
    result = json.loads((tmp_path/'completion.json').read_text())
    assert result['returncode'] == 0
    assert result['status'] == 'namespace_process_returned_physical_review_required'
