"""A live-worker resource limit must never signal a stale or unrelated owner."""
import importlib.util
import json
import os
from pathlib import Path
import signal
import socket
import subprocess
import sys
from types import SimpleNamespace

import psutil
import pytest

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location('owned_host_memory', ROOT/'scripts/fleet/limit_owned_host_memory.py')
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def identity(process):
    return dict(host=socket.gethostname(), pid=process.pid,
                process_start=Path(f'/proc/{process.pid}/stat').read_text().rsplit(')', 1)[1].split()[19])


def test_stale_or_mismatched_owner_cannot_be_signaled(tmp_path):
    campaign = tmp_path/'case'
    program = tmp_path/'worker.py'
    program.write_text('import time; time.sleep(30)\n')
    process = subprocess.Popen([sys.executable, str(program), str(campaign)])
    try:
        owner = identity(process)
        for bad in (dict(owner, process_start='invalid'), dict(owner, host='another-host')):
            with pytest.raises(ValueError):
                MODULE.enforce(bad, campaign, program, 1, 1)
            assert process.poll() is None
        with pytest.raises(ValueError):
            MODULE.enforce(owner, campaign, tmp_path/'other.py', 1, 1)
        assert process.poll() is None
    finally:
        process.terminate()
        process.wait(timeout=10)


def test_memory_stop_targets_exact_worker_and_preserves_checkpoint(tmp_path, monkeypatch):
    program, campaign = tmp_path/'worker.py', tmp_path/'case'
    campaign.mkdir()
    program.write_text('''import signal,time,sys
from pathlib import Path
case=Path(sys.argv[1])
def stop(signum,frame):
    (case/'checkpoint.json').write_text('{"interrupted":true,"qualified":false}')
    raise SystemExit(0)
signal.signal(signal.SIGTERM,stop)
(case/'ready').touch()
while True:time.sleep(.02)
''')
    process = subprocess.Popen([sys.executable, str(program), str(campaign)])
    try:
        import time
        deadline = time.monotonic()+10
        while not (campaign/'ready').exists():
            assert process.poll() is None and time.monotonic() < deadline
            time.sleep(.01)
        owner = identity(process)
        state = MODULE.enforce(owner, campaign, program, 2**50, 1)
        assert state['status'] == 'owned_worker_within_memory_bound'
        assert process.poll() is None
        state = MODULE.enforce(owner, campaign, program, 1, 1)
        assert state['status'] == 'owned_worker_memory_stop_requested'
        process.wait(timeout=10)
        assert json.loads((campaign/'checkpoint.json').read_text()) == dict(interrupted=True, qualified=False)
    finally:
        if process.poll() is None:
            process.terminate()
            process.wait(timeout=10)
