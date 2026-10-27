import importlib.util
import os
from pathlib import Path
import pytest

ROOT=Path(__file__).resolve().parents[1]
s=importlib.util.spec_from_file_location('observe_polish',ROOT/'scripts/pipeline/observe_counterpoise_polish.py');m=importlib.util.module_from_spec(s);s.loader.exec_module(m)


def test_process_start_identity_is_required_for_a_live_owner():
    start=Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()[19]
    owner=dict(host=os.uname().nodename,pid=os.getpid(),process_start=start)
    assert m.alive(owner)
    assert not m.alive(dict(owner,process_start='wrong-instance'))
    assert not m.alive(dict(owner,pid=999999999))


def test_unverified_remote_owner_is_not_assumed_stopped():
    with pytest.raises(ValueError):m.alive(dict(host='unknown-host',pid=1,process_start='unknown'))
