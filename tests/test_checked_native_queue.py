import importlib.util
import json
from pathlib import Path
import subprocess

import pytest

ROOT = Path(__file__).resolve().parents[1]
definition = importlib.util.spec_from_file_location('checked_native_queue',
    ROOT/'scripts/pipeline/run_checked_native_queue.py')
queue = importlib.util.module_from_spec(definition)
definition.loader.exec_module(queue)


def preparation(tmp_path):
    case = tmp_path/'case'
    case.mkdir()
    (case/'plan.json').write_text('original settings')
    (case/'input.xyz').write_text('original coordinates')
    spec = dict(key='task', attempt=str(case), maximum_host_RSS_GiB=80,
        plan_sha256=queue.sha256_file(case/'plan.json'), geometry_sha256=queue.sha256_file(case/'input.xyz'))
    namespace = tmp_path/'namespace.json'
    namespace.write_text(json.dumps(spec))
    row = dict(status='native_inputs_prepared', key='task', energy_namespace='declared',
        native_namespace=str(namespace), native_namespace_sha256=queue.sha256_file(namespace))
    return dict(entries=[row]), spec


def test_exact_growing_preparation_stream_retains_original_settings(tmp_path):
    stream, spec = preparation(tmp_path)
    bindings = {}
    assert queue.ready_inputs(stream, {'declared'}, bindings)[0][1] == spec
    assert bindings == {stream['entries'][0]['native_namespace']: stream['entries'][0]['native_namespace_sha256']}


@pytest.mark.parametrize('mutation', ['undeclared', 'duplicate', 'namespace', 'geometry', 'prior_binding'])
def test_changed_or_duplicate_input_cannot_be_dispatched(tmp_path, mutation):
    stream, spec = preparation(tmp_path)
    row = stream['entries'][0]
    bindings = {}
    if mutation == 'undeclared':row['energy_namespace'] = 'different'
    if mutation == 'duplicate':stream['entries'].append(dict(row))
    if mutation == 'namespace':Path(row['native_namespace']).write_text('{}')
    if mutation == 'geometry':(Path(spec['attempt'])/'input.xyz').write_text('wrong geometry')
    if mutation == 'prior_binding':bindings[row['native_namespace']] = 'different'
    with pytest.raises(ValueError):
        queue.ready_inputs(stream, {'declared'}, bindings)


@pytest.mark.parametrize('liveness,expected', [(True, 'running'), (False, 'stopped_attempt_preserved'),
                                            (None, 'ownership_unverified_preserved')])
def test_live_dead_and_unknown_owners_never_reopen_scientific_outputs(tmp_path, liveness, expected):
    stage = tmp_path/'experimental_native'
    stage.mkdir()
    (stage/'status.json').write_text(json.dumps(dict(status='running', owner=dict(pid=123))))
    def alive(owner):
        if liveness is None:raise subprocess.TimeoutExpired('ssh', 20)
        return liveness
    assert queue.inspect(dict(attempt=str(tmp_path)), alive)[0] == expected


def test_guard_and_submission_claim_block_duplicate_starts(tmp_path):
    claim = tmp_path/'native_dispatch_claim'
    claim.mkdir()
    (claim/'owner.json').write_text('{"host":"node"}')
    assert queue.inspect(dict(attempt=str(tmp_path)), lambda _: False)[0] == 'submission_claim_preserved'
    guard = tmp_path/'residency_guard'
    guard.mkdir()
    (guard/'owner.json').write_text('{"host":"node"}')
    assert queue.inspect(dict(attempt=str(tmp_path)), lambda _: False)[0] == 'guard_present_start_pending'


def test_finished_worker_still_needs_closed_ownership_before_review(tmp_path):
    stage = tmp_path/'experimental_native'
    stage.mkdir()
    (stage/'status.json').write_text(json.dumps(dict(owner=dict(pid=1))))
    spec = dict(attempt=str(tmp_path))
    assert queue.physical_review(tmp_path/'namespace.json', spec, tmp_path, lambda _: True) is None
    assert queue.physical_review(tmp_path/'namespace.json', spec, tmp_path, lambda _: False) is None


def gate_refusal(tmp_path):
    case, raw, archive, output, placement = [tmp_path/name for name in ('case', 'raw', 'archive', 'queue', 'placement')]
    for p in (case, raw, archive, output, placement):p.mkdir()
    spec = dict(key='key', attempt=str(case), raw_parent=str(raw), archive=str(archive), prior_root='/old',
        original_campaign='/old/results/polymer_models/campaigns/native/key',
        original_raw_parent='/old/calcs/polymer_models/native', original_raw_root='/old/calcs/polymer_models/native/key')
    claim = case/'native_dispatch_claim'
    claim.mkdir()
    (claim/'owner.json').write_text(json.dumps(dict(queue=str(output))))
    (claim/'launch.json').write_text(json.dumps(dict(status='controller_spawned_scientific_start_unverified',
        owner=dict(host='node', pid=1, process_start='first'), argv=['program', '--out', str(placement)])))
    (placement/'completion.json').write_text(json.dumps(dict(returncode=1, GPU=0, host='node', pid=2, process_start='second')))
    (placement/'command.log').write_text('gpu_reserve: gpu0 BUSY/reserved -> not running\n')
    return spec, output, claim


def test_proven_gate_refusal_retries_without_consuming_a_scientific_attempt(tmp_path):
    spec, output, claim = gate_refusal(tmp_path)
    assert queue.release_gate_refusal(spec, output, lambda _: False)
    assert not claim.exists()
    proof = list((output/'gate_contention_history/key').glob('*/resolution.json'))
    assert len(proof) == 1 and json.loads(proof[0].read_text())['inference_attempt_consumed'] is False


@pytest.mark.parametrize('uncertainty', ['live', 'unreachable', 'scientific_output'])
def test_gate_refusal_never_reopens_an_uncertain_or_started_worker(tmp_path, uncertainty):
    spec, output, claim = gate_refusal(tmp_path)
    def alive(owner):
        if uncertainty == 'unreachable':raise subprocess.TimeoutExpired('ssh', 20)
        return uncertainty == 'live'
    if uncertainty == 'scientific_output':(Path(spec['attempt'])/'experimental_native').mkdir()
    assert not queue.release_gate_refusal(spec, output, alive)
    assert claim.exists()
