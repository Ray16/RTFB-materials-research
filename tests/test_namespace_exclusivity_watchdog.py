"""An actual owned TERM-ignoring subprocess is stopped; wrong identity is safe."""
import importlib.util
from pathlib import Path
import subprocess
import sys

import pytest
from redox.core.common import ROOT


def watcher():
    path=ROOT/'scripts/fleet/watch_namespace_gpu_owners.py'
    spec=importlib.util.spec_from_file_location('test_namespace_watchdog',path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module


def process(plan):
    code='import signal,time;signal.signal(signal.SIGTERM,signal.SIG_IGN);print("ready",flush=True);time.sleep(60)'
    child=subprocess.Popen([sys.executable,'-c',code,str(plan)],stdout=subprocess.PIPE,text=True)
    assert child.stdout.readline().strip()=='ready'
    return child


def test_owned_worker_ignoring_term_is_escalated_without_other_process_signals(tmp_path):
    monitor=watcher();plan=tmp_path/'plan.json';child=process(plan)
    try:
        result=monitor.stop_owned(child.pid,monitor.process_start(child.pid),plan,grace=.1)
        assert result==dict(stopped=True,forced=True)
        # psutil.wait reaps this test's child before Popen can recover its exit
        # status; the authoritative evidence is that the exact PID is gone.
        assert not Path(f'/proc/{child.pid}').exists()
        child.wait(timeout=2)
    finally:
        if child.poll() is None:child.kill();child.wait()


def test_wrong_start_or_command_never_signals_a_process(tmp_path):
    monitor=watcher();plan=tmp_path/'plan.json';child=process(plan)
    try:
        with pytest.raises(ValueError):monitor.stop_owned(child.pid,'wrong-start',plan,grace=.1)
        with pytest.raises(ValueError):monitor.stop_owned(child.pid,monitor.process_start(child.pid),tmp_path/'different.json',grace=.1)
        assert child.poll() is None
    finally:child.kill();child.wait()
