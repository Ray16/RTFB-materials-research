"""Residency guard stops only the exact owned worker and never a peer process."""
import importlib.util
from pathlib import Path
import subprocess
import sys

import pytest

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('GPU_exclusivity_watch',ROOT/'scripts/fleet/watch_gpu_exclusivity.py')
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def test_any_foreign_resident_on_assigned_UUID_is_a_conflict():
    assert m.conflicts(10,3,{3:'ours',4:'elsewhere'},[(10,'ours'),(11,'ours'),(12,'elsewhere')])==[11]
    assert m.conflicts(10,3,{3:'ours'},[(10,'ours'),(10,'ours')])==[]
    with pytest.raises(ValueError):
        m.conflicts(10,3,{},[])


def test_only_exact_owned_process_is_stopped_and_checkpoint_is_preserved(tmp_path):
    checkpoint=tmp_path/'saved.xyz';checkpoint.write_text('last complete evaluated geometry')
    command=[sys.executable,'-c','import time;time.sleep(30)',str(tmp_path)]
    ours=subprocess.Popen(command)
    peer=subprocess.Popen(command)
    try:
        start=m.process_start(ours.pid)
        with pytest.raises(ValueError):
            m.stop_owned(ours.pid,start+'0',tmp_path)
        assert ours.poll() is None and peer.poll() is None
        with pytest.raises(ValueError):
            m.stop_owned(ours.pid,start,tmp_path/'different_campaign')
        assert ours.poll() is None
        result=m.stop_owned(ours.pid,start,tmp_path)
        assert result['stopped'] and ours.poll() is not None
        assert peer.poll() is None
        assert checkpoint.read_text()=='last complete evaluated geometry'
    finally:
        for p in (ours,peer):
            if p.poll() is None:p.terminate()
            p.wait(timeout=5)
