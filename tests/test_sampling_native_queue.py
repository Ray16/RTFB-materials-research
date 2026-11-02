"""A native queue keeps the full objective and cannot consume a partial stream."""
import json
from pathlib import Path
import pytest

from redox.core.hashing import sha256_file
from redox.workflows import sampling_native_queue as api


def stream():
    rows = [dict(key=f'key{i}', status='pending') for i in range(240)]
    manifest = dict(target_states=240, experimental_method_id='method', states=rows)
    result = dict(status='complete_roster_corrected_electronic_review_snapshot', declared_inputs=240,
                  experimental_method_id='method', independent_review_errors=[], states=[dict(r) for r in rows])
    return manifest, result


@pytest.mark.parametrize('bad', ['missing', 'duplicate', 'extra', 'method', 'errors'])
def test_incomplete_or_different_electronic_stream_cannot_release_native_inputs(bad):
    manifest, result = stream()
    if bad == 'missing': result['states'].pop()
    if bad == 'duplicate': result['states'][-1] = dict(result['states'][0])
    if bad == 'extra': result['states'][-1]['key'] = 'undeclared'
    if bad == 'method': result['experimental_method_id'] = 'other'
    if bad == 'errors': result['independent_review_errors'] = [dict(key='key1')]
    with pytest.raises(ValueError): api.checked_stream(manifest, result)


def test_pending_inputs_are_retained_without_receiving_native_values(tmp_path):
    manifest, result = stream()
    result['states'][0]['status'] = 'electronic_independently_checked'
    assert len(api.checked_stream(manifest, result)) == 1
    path = tmp_path / 'manifest.json'
    path.write_text(json.dumps(manifest))
    collected = api.collect(path, [], tmp_path / 'reviews')
    assert collected['declared_inputs'] == len(collected['states']) == 240
    assert collected['checked_native_inputs'] == 0
    assert all(r['G_by_cutoff_eV'] is None and not r['native_G_qualified'] for r in collected['states'])
    assert not collected['occupied_ensemble_cycles_qualified']


def test_preparation_requires_largest_pilot_and_source_bound_energy(tmp_path, monkeypatch):
    manifest, result = stream()
    manifest_path = tmp_path / 'manifest.json'
    manifest_path.write_text(json.dumps(manifest))
    result['states'][0].update(status='electronic_independently_checked', review=str(tmp_path/'missing'), review_sha256='bad')
    with pytest.raises(ValueError, match='largest complete native pilot'):
        api.prepare_ready(manifest_path, result, [], None, None, None, tmp_path/'native')


def test_peer_extensions_preserve_existing_namespaces_and_reject_mutations(tmp_path):
    first, second = tmp_path/'one.json', tmp_path/'two.json'
    first.write_text('{}'); second.write_text('{}')
    path = tmp_path/'peers.json'
    initial = dict(namespace=str(first), namespace_sha256=sha256_file(first))
    path.write_text(json.dumps([initial]))
    row = dict(namespace=str(second), namespace_sha256=sha256_file(second))
    assert api.extend_peers(path, [row]) == [initial, row]
    assert api.extend_peers(path, [row]) == [initial, row]
    second.write_text('changed')
    with pytest.raises(ValueError): api.extend_peers(path, [row])


@pytest.mark.parametrize('condition,tries,terminal', [
    ('scientific_failure_preserved', 1, True),
    ('stopped_attempt_preserved', 1, True),
    ('prelaunch_refusal_no_scientific_attempt', 3, True),
    ('prelaunch_refusal_no_scientific_attempt', 2, False),
    ('running', 1, False), ('ownership_unverified', 1, False),
    ('completed_review_pending', 1, False)])
def test_terminal_failures_do_not_masquerade_as_a_live_native_queue(condition, tries, terminal):
    import importlib.util
    from redox.core.common import ROOT
    definition = importlib.util.spec_from_file_location('native_queue_terminal_test', ROOT/'scripts/fleet/run_corrected_native_queue.py')
    queue = importlib.util.module_from_spec(definition)
    definition.loader.exec_module(queue)
    rows = [dict(key='one', status=condition, submission_receipts=tries)]
    assert queue.terminal_incomplete(240, rows, [], set()) is terminal
    assert not queue.terminal_incomplete(239, rows, [], set())
    assert not queue.terminal_incomplete(240, rows, ['ready'], set())
