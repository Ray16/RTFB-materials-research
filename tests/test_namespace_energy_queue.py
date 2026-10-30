import importlib.util
import json
from pathlib import Path
import subprocess

import pytest

ROOT = Path(__file__).resolve().parents[1]
specification = importlib.util.spec_from_file_location('namespace_energy_queue_contract',
    ROOT/'scripts/pipeline/run_namespace_energy_queue.py')
queue = importlib.util.module_from_spec(specification)
specification.loader.exec_module(queue)


def resources(**updates):
    result = dict(GPU_MiB=32768, available=400*2**30, promised=100*2**30, load=2, cpus=80)
    result.update(updates)
    return result


def test_host_memory_promises_count_even_when_observed_rss_is_small():
    state = dict(maximum_host_RSS_GiB=160)
    assert queue.fits(state, resources())
    assert not queue.fits(state, resources(promised=250*2**30))


@pytest.mark.parametrize('updates', [dict(GPU_MiB=16384), dict(load=79), dict(available=128*2**30)])
def test_unsuitable_capacity_never_enters_submission(updates):
    assert not queue.fits(dict(maximum_host_RSS_GiB=160), resources(**updates))


def test_guard_before_scientific_status_prevents_duplicate_launch(tmp_path):
    guard = tmp_path/'residency_guard'
    guard.mkdir()
    (guard/'owner.json').write_text(json.dumps(dict(host='worker', pid=42, process_start='start')))
    state, _ = queue.inspect_state(dict(attempt=str(tmp_path)), lambda owner: False)
    assert state == 'guard_present_scientific_start_pending'


def test_unreachable_scientific_owner_is_preserved(tmp_path):
    stage = tmp_path/'experimental_energy'
    stage.mkdir()
    (stage/'status.json').write_text(json.dumps(dict(status='running', owner=dict(host='worker', pid=42))))

    def uncertain(owner):
        raise subprocess.TimeoutExpired('ssh', 20)

    state, _ = queue.inspect_state(dict(attempt=str(tmp_path)), uncertain)
    assert state == 'ownership_unverified_preserved'


def test_dead_owner_does_not_reopen_partial_scientific_outputs(tmp_path):
    stage = tmp_path/'experimental_energy'
    stage.mkdir()
    (stage/'status.json').write_text(json.dumps(dict(status='running', owner=dict(host='worker', pid=42))))
    state, _ = queue.inspect_state(dict(attempt=str(tmp_path)), lambda owner: False)
    assert state == 'stopped_attempt_preserved'
