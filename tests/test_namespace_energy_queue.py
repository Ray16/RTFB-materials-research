import importlib.util
import json
from pathlib import Path
import subprocess
from types import SimpleNamespace

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


def test_node_timeout_does_not_discard_verified_peer_availability(monkeypatch):
    def query(argv, **kwargs):
        if argv[-1] == 'unreachable':
            raise subprocess.TimeoutExpired(argv, 15)
        return SimpleNamespace(stdout='healthy 2\n')

    monkeypatch.setattr(queue.subprocess, 'run', query)
    slots, errors = queue.discover_slots(['healthy', 'unreachable'])
    assert slots == [('healthy', 2)]
    assert errors[0]['host'] == 'unreachable'
    assert errors[0]['availability'] == 'unverified_no_submission'


def test_malformed_availability_cannot_authorize_another_host(monkeypatch):
    monkeypatch.setattr(queue.subprocess, 'run', lambda *a, **k: SimpleNamespace(stdout='different-host 0\n'))
    slots, errors = queue.discover_slots(['requested-host'])
    assert not slots and len(errors) == 1


def refused_claim(tmp_path):
    case, raw, archive, output, placement = [tmp_path/name for name in ('case', 'raw', 'archive', 'queue', 'placement')]
    for path in (case, raw, archive, output, placement):
        path.mkdir()
    spec = dict(key='key', attempt=str(case), raw_parent=str(raw), archive=str(archive), prior_root='/old',
                original_campaign='/old/results/polymer_models/campaigns/batch/key',
                original_raw_parent='/old/calcs/polymer_models/batch',
                original_raw_root='/old/calcs/polymer_models/batch/key')
    claim = case/'dispatch_claim'
    claim.mkdir()
    (claim/'owner.json').write_text(json.dumps(dict(queue=str(output))))
    receipt = dict(status='controller_spawned_scientific_start_unverified',
                   owner=dict(host='node', pid=101, process_start='first'), argv=['program', '--out', str(placement)])
    (claim/'launch.json').write_text(json.dumps(dict(receipts=[receipt])))
    (placement/'completion.json').write_text(json.dumps(dict(returncode=1, GPU=0, host='node', pid=102, process_start='second')))
    (placement/'command.log').write_text('gpu_reserve: gpu0 BUSY/reserved -> not running\n')
    return spec, output, claim


def test_proven_stopped_gate_refusal_is_contention_not_a_scientific_attempt(tmp_path):
    spec, output, claim = refused_claim(tmp_path)
    assert queue.release_gate_contention(spec, output, lambda owner: False)
    assert not claim.exists()
    proof = list((output/'gate_contention_history/key').glob('*/resolution.json'))
    assert len(proof) == 1 and json.loads(proof[0].read_text())['inference_attempt_consumed'] is False


def test_live_owner_or_scientific_output_blocks_contention_release(tmp_path):
    spec, output, claim = refused_claim(tmp_path)
    assert not queue.release_gate_contention(spec, output, lambda owner: True)
    (Path(spec['raw_parent'])/'key').mkdir()
    assert not queue.release_gate_contention(spec, output, lambda owner: False)
    assert claim.exists()


def test_unknown_owner_blocks_contention_release(tmp_path):
    spec, output, claim = refused_claim(tmp_path)

    def uncertain(owner):
        raise subprocess.TimeoutExpired('ssh', 20)

    assert not queue.release_gate_contention(spec, output, uncertain)
    assert claim.exists()
