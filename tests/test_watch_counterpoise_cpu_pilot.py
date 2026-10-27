import importlib.util
from pathlib import Path

import pytest


ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('cpu_pilot_watcher',ROOT/'scripts/pipeline/watch_counterpoise_cpu_pilot.py')
watcher=importlib.util.module_from_spec(spec);spec.loader.exec_module(watcher)


def test_review_waits_for_converged_status_and_owner_exit():
    running={'status':'CPU_pilot_SCF_running_unqualified'}
    completed={'status':'CPU_pilot_converged_independent_review_pending'}
    assert watcher.classify(running,'live')=='waiting_for_pilot'
    assert watcher.classify(running,'stopped')=='pilot_owner_stopped_without_terminal_status'
    assert watcher.classify(completed,'live')=='waiting_for_owner_exit'
    assert watcher.classify(completed,'stopped')=='review_ready'


def test_failed_or_unknown_pilot_is_not_promoted():
    assert watcher.classify({'status':'CPU_pilot_interrupted_or_failed_preserved'},'stopped')=='pilot_failed_preserved'
    with pytest.raises(ValueError,match='Unknown pilot state'):
        watcher.classify({'status':'invented'},'stopped')
