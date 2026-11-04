import importlib.util
import json
import os
from pathlib import Path
import sys

import pytest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('CP_checkpoint_diagnostic',
    ROOT/'scripts/validation/diagnose_archived_counterpoise_checkpoint.py')
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def argv(tmp_path, extra=()):
    archive = tmp_path/'archive'
    archive.mkdir()
    args = ['CP-diagnostic', '--archive', str(archive), '--prior-root', '/old',
            '--checkpoint-root', '/old/calcs/checkpoint', '--out', str(tmp_path/'out'), *extra]
    return archive, args


def test_existing_diagnostic_output_is_never_overwritten(tmp_path, monkeypatch):
    _, args = argv(tmp_path)
    (tmp_path/'out').mkdir()
    monkeypatch.setattr(sys, 'argv', args)
    with pytest.raises(FileExistsError):
        module.main()


@pytest.mark.parametrize('extra', [('--threads', '0'), ('--maximum-rss-gib', '64')])
def test_invalid_or_changed_resource_request_stops_before_any_computation(tmp_path, monkeypatch, extra):
    _, args = argv(tmp_path, extra)
    monkeypatch.setattr(sys, 'argv', args)
    with pytest.raises(ValueError, match='retained 32 GiB'):
        module.main()
    assert not (tmp_path/'out').exists()


def test_live_checkpoint_producer_cannot_be_read_as_a_stopped_attempt(tmp_path, monkeypatch):
    archive, args = argv(tmp_path)
    raw = archive/'calcs/checkpoint'
    raw.mkdir(parents=True)
    (raw/'launch.json').write_text(json.dumps(dict(role='AB')))
    owner = dict(host=os.uname().nodename, pid=os.getpid(),
                 process_start=Path('/proc/self/stat').read_text().rsplit(')', 1)[1].split()[19])
    (raw/'status.json').write_text(json.dumps(dict(role='AB', owner=owner)))
    monkeypatch.setattr(sys, 'argv', args)
    with pytest.raises(ValueError, match='still live'):
        module.main()
    assert not (tmp_path/'out').exists()
