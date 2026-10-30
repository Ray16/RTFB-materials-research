import importlib.util
import json
import os
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
definition = importlib.util.spec_from_file_location('energy_native_priority',
    ROOT/'scripts/pipeline/run_energy_queue_with_native_priority.py')
priority = importlib.util.module_from_spec(definition)
definition.loader.exec_module(priority)


def test_live_exact_dispatcher_with_pending_states_receives_priority(tmp_path):
    path = tmp_path/'status.json'
    record = dict(owner=dict(host=os.uname().nodename, pid=os.getpid(),
        process_start=Path('/proc/self/stat').read_text().rsplit(')', 1)[1].split()[19]), pending=2)
    path.write_text(json.dumps(record))
    assert priority.native_waiting(path)
    record['pending'] = 0
    path.write_text(json.dumps(record))
    assert not priority.native_waiting(path)


@pytest.mark.parametrize('change', ['pid_reuse', 'terminal', 'missing_process'])
def test_stopped_or_reused_owner_cannot_starve_the_energy_queue(tmp_path, change):
    path = tmp_path/'status.json'
    record = dict(owner=dict(host=os.uname().nodename, pid=os.getpid(),
        process_start=Path('/proc/self/stat').read_text().rsplit(')', 1)[1].split()[19]), pending=2)
    if change == 'pid_reuse':record['owner']['process_start'] = 'different'
    if change == 'terminal':record['finished'] = 1
    if change == 'missing_process':record['owner']['pid'] = 2**30
    path.write_text(json.dumps(record))
    assert not priority.native_waiting(path)


def test_missing_priority_state_is_not_silently_assumed_idle(tmp_path):
    with pytest.raises(ValueError):
        priority.native_waiting(tmp_path/'missing.json')
