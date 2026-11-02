"""Native values cannot survive a changed component, state or thermal context."""
import json
from pathlib import Path
import pytest

from redox.core.hashing import sha256_file
from redox.validation import sampling_native_receipts as api


def fixture(tmp_path, monkeypatch):
    case = tmp_path / 'case'
    raw = case / 'raw'
    raw.mkdir(parents=True)
    namespace = tmp_path / 'namespace.json'
    namespace.write_text(json.dumps(dict(case=str(case))))
    geometry = case / 'input.xyz'
    geometry.write_text('exact geometry test fixture')
    task = dict(key='key', charge=0, mult=1)
    model, protocol = dict(name='test'), dict(xc='test')
    thermal = dict(pressure='specified translation convention', temperature_K=298.15)
    electronic = dict(namespace=str(namespace), task=task, E_smd_eV=-20.,
                      experimental_model=model, electronic_protocol=protocol, source_sha256={})
    electronic_path = case / 'electronic.json'
    electronic_path.write_text(json.dumps(electronic))
    partition = dict(nuclear=dict(nuclear_by_cutoff_eV={k: .1 for k in ('25', '50', '100')}),
                     G_by_cutoff_eV={k: -19.9 for k in ('25', '50', '100')}, spin_partition_eV=0.,
                     experimental_method_id='test', derivative_checks_passed=True, primary_candidate_admission=False)
    partition_path = raw / 'experimental_partition.json'
    partition_path.write_text(json.dumps(partition))
    plan = dict(role='corrected_sampling_native', raw_root=str(raw), electronic_review=str(electronic_path),
                tasks=[task], experimental_model=model, electronic_protocol=protocol,
                experimental_native_settings=dict(xc='native'), thermal_protocol=thermal,
                scored_geometry_sha256=sha256_file(geometry), source_sha256={})
    plan_path = case / 'plan.json'
    plan_path.write_text(json.dumps(plan))
    report = dict(status='corrected_sampling_native_independently_reconstructed', namespace=str(namespace),
                  key='key', experimental_method_id='test', native_G_qualified=True,
                  primary_candidate_admission=False, conformer_ensemble_qualified=False,
                  task=task, charge=0, multiplicity=1, experimental_model=model, electronic_protocol=protocol,
                  native_settings=plan['experimental_native_settings'], thermal_protocol=thermal,
                  temperature_K=298.15, standard_state=thermal['pressure'], geometry_sha256=sha256_file(geometry),
                  thermal_source=partition['nuclear'], electronic_E_smd_eV=-20., spin_partition_eV=0.,
                  G_by_cutoff_eV=partition['G_by_cutoff_eV'], nuclear_by_cutoff_eV=partition['nuclear']['nuclear_by_cutoff_eV'],
                  minimum_frequency_cm=10., maximum_column_error=.0001, observed_peak_GPU_MiB=20000,
                  peak_host_RSS_bytes=30*2**30, actual_SCF_frames_checked=17,
                  source_sha256={str(p): sha256_file(p) for p in (namespace, plan_path, partition_path, electronic_path)})
    receipt = tmp_path / 'review.json'
    receipt.write_text(json.dumps(report))
    monkeypatch.setattr(api, 'verify', lambda p: dict(case=str(case)))
    monkeypatch.setattr(api, 'checked_energy', lambda *args: electronic)
    return receipt, namespace, report, plan, partition


def test_native_receipt_preserves_own_components_and_all_cutoffs(tmp_path, monkeypatch):
    path, namespace, report, _, _ = fixture(tmp_path, monkeypatch)
    assert api.checked_native(path, namespace, 'key', 'test')['G_by_cutoff_eV'] == report['G_by_cutoff_eV']


@pytest.mark.parametrize('bad', ['G', 'nuclear', 'spin', 'charge', 'multiplicity', 'thermal', 'standard_state',
                               'temperature', 'model', 'geometry', 'cutoff', 'missing_raw', 'ensemble', 'method'])
def test_native_receipt_refuses_context_or_value_substitution(tmp_path, monkeypatch, bad):
    path, namespace, report, _, _ = fixture(tmp_path, monkeypatch)
    if bad == 'G': report['G_by_cutoff_eV']['50'] += .01
    if bad == 'nuclear': report['nuclear_by_cutoff_eV']['50'] += .01
    if bad == 'spin': report['spin_partition_eV'] = .1
    if bad == 'charge': report['charge'] = 1
    if bad == 'multiplicity': report['multiplicity'] = 2
    if bad == 'thermal': report['thermal_protocol']['pressure'] = 'other'
    if bad == 'standard_state': report['standard_state'] = 'other'
    if bad == 'temperature': report['temperature_K'] = 300
    if bad == 'model': report['experimental_model'] = dict(name='other')
    if bad == 'geometry': (tmp_path / 'case/input.xyz').write_text('different')
    if bad == 'cutoff': report['G_by_cutoff_eV'].pop('25')
    if bad == 'missing_raw': report['source_sha256'].pop(str(tmp_path / 'case/raw/experimental_partition.json'))
    if bad == 'ensemble': report['conformer_ensemble_qualified'] = True
    if bad == 'method': report['experimental_method_id'] = 'other'
    path.write_text(json.dumps(report))
    with pytest.raises(ValueError): api.checked_native(path, namespace, 'key', 'test')


def test_raw_and_report_can_agree_but_wrong_electronic_nuclear_sum_is_rejected(tmp_path, monkeypatch):
    path, namespace, report, _, partition = fixture(tmp_path, monkeypatch)
    partition['G_by_cutoff_eV']['50'] += 1
    report['G_by_cutoff_eV']['50'] = partition['G_by_cutoff_eV']['50']
    raw = tmp_path / 'case/raw/experimental_partition.json'
    raw.write_text(json.dumps(partition))
    report['source_sha256'][str(raw)] = sha256_file(raw)
    path.write_text(json.dumps(report))
    with pytest.raises(ValueError, match='sum'): api.checked_native(path, namespace, 'key', 'test')


def test_cache_invalidates_same_length_backing_changes(tmp_path, monkeypatch):
    path, namespace, report, _, _ = fixture(tmp_path, monkeypatch)
    cache = api.ReceiptCache()
    # The fixture uses one namespace for both components. Its plan still has
    # the complete binding fields needed to exercise stat invalidation.
    calls = []
    original = api.checked_native
    def counted(*args):
        calls.append(1)
        return original(*args)
    monkeypatch.setattr(api, 'checked_native', counted)
    cache.read(path, namespace, 'key', 'test')
    cached = cache.read(path, namespace, 'key', 'test')
    cached['G_by_cutoff_eV']['50'] = 5
    assert len(calls) == 1
    assert cache.read(path, namespace, 'key', 'test')['G_by_cutoff_eV']['50'] == -19.9
    raw = tmp_path / 'case/raw/experimental_partition.json'
    text = raw.read_text()
    raw.write_text(text.replace('-19.9', '-18.9'))
    with pytest.raises(ValueError): cache.read(path, namespace, 'key', 'test')
    assert len(calls) == 2


def test_cached_receipt_detects_a_later_residency_interruption(tmp_path, monkeypatch):
    path, namespace, _, _, _ = fixture(tmp_path, monkeypatch)
    cache = api.ReceiptCache()
    cache.read(path, namespace, 'key', 'test')
    marker = tmp_path/'case/residency_guard/external_interruption.json'
    marker.parent.mkdir()
    marker.write_text('{}')
    with pytest.raises(ValueError, match='residency'): cache.read(path, namespace, 'key', 'test')
