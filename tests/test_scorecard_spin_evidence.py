"""A helper energy alone cannot qualify a selected spin for a primary ranking."""
import json

import pytest

from redox.core import common, protocol
from redox.qm import dft
from redox.screening import scorecard
from redox.validation import spins


@pytest.fixture
def helper(tmp_path, monkeypatch):
    directory = tmp_path / 'spincheck/m/red2_m3'
    directory.mkdir(parents=True)
    (directory / 'opt.xyz').write_text('2\n\nH 0 0 0\nH 0 0 .8\n')
    raw = dict(charge=0, mult=3, g_thermal_eV=0, n_imag=0)
    (directory / 'result.json').write_text(json.dumps(raw))
    address = protocol.record_path(directory, 0, 3)
    address.parent.mkdir()
    address.write_text(json.dumps(dict(protocol=protocol.ACTIVE_SP, protocol_hash=address.stem,
        geom_sha1=protocol.geom_sha1(directory / 'opt.xyz'), charge=0, mult=3,
        converged_gas=True, converged_smd=True, e_gas_eV=-8, e_smd_eV=-9)))
    identity = {'geometry': 'current'}
    report = dict(status='spin_gap_verified', chosen_identity=identity.copy(),
                  alternative_mult=3, electronic_gap_eV=1)
    monkeypatch.setattr(scorecard, 'SPINCHECK', tmp_path / 'spincheck')
    monkeypatch.setattr(scorecard, 'UMA', tmp_path / 'uma')
    monkeypatch.setattr(common, 'read_result', lambda *args, **kw:
        dict(charge=0, mult=1, e_smd_eV=-10, g_thermal_eV=0, n_imag=0))
    monkeypatch.setattr(spins, 'snapshot', lambda *args: identity)
    monkeypatch.setattr(spins, 'read_comparison', lambda *args: report)
    monkeypatch.setattr(dft, '_opt_protocol_mismatch', lambda *args: '')
    return directory, raw, address, identity, report


@pytest.mark.parametrize('mutation', [
    {'status': 'running'}, {'status': 'failed'}, {'chosen_identity': {'geometry': 'old'}},
    {'electronic_gap_eV': .8}, {'electronic_gap_eV': float('nan')},
    {'alternative_mult': 1}, {'ensemble_primary_integration_pending': True},
])
def test_unfinished_or_stale_helper_never_overrides_spin_confidence(helper, mutation):
    helper[-1].update(mutation)
    gap, source = scorecard._min_spin_gap_eV('m')
    assert gap is None and source.startswith('pending:')


def test_unqualified_helper_status_blocks_an_otherwise_completed_report(helper):
    directory, raw, *_ = helper
    raw['spin_validation_status'] = 'failed'
    (directory / 'result.json').write_text(json.dumps(raw))
    assert scorecard._min_spin_gap_eV('m') == (None, 'pending:alternate_spin_validation_pending:red2')


def test_close_electronic_gap_requires_matching_native_thermochemistry(helper):
    _, _, address, _, report = helper
    record = json.loads(address.read_text())
    record['e_smd_eV'] = -9.85
    address.write_text(json.dumps(record))
    report['electronic_gap_eV'] = .15
    assert scorecard._min_spin_gap_eV('m') == (None, 'pending:matching_native_spin_thermals_required:red2')


def test_completed_isolated_comparison_retains_explicit_thermal_uncertainty(helper):
    assert scorecard._min_spin_gap_eV('m') == (1, 'dft_electronic(thermal_unvalidated):red2')


def test_live_spin_writer_blocks_ranking_even_during_helper_directory_swap(tmp_path, monkeypatch):
    import fcntl
    from redox.validation import spin_evidence
    monkeypatch.setattr(common, 'ROOT', tmp_path)
    locks=tmp_path/'logs/pipeline_locks'
    locks.mkdir(parents=True)
    with (locks/'m_spincheck.lock').open('w') as owner:
        fcntl.flock(owner,fcntl.LOCK_EX|fcntl.LOCK_NB)
        assert spin_evidence.spin_comparison_running('m')
        assert scorecard._min_spin_gap_eV('m')==(None,'pending:active_spin_comparison')
        assert not spin_evidence.spin_comparison_running('other')
    assert not spin_evidence.spin_comparison_running('m')  # an old lock file is not a live owner
def test_native_spin_gap_requires_current_corrected_partition_files(tmp_path):
    import hashlib
    from redox.validation.spin_evidence import native_partition_dependencies_ok
    chosen,alternative=tmp_path/'chosen.json',tmp_path/'alternative.json'
    chosen.write_text('{}');alternative.write_text('{"spin":3}')
    protocol={'schema':1,'model':'qualified projected partition'}
    left=dict(thermal_partition_protocol=protocol,thermal_postprocessing_path=str(chosen))
    right=dict(thermal_partition_protocol=protocol,thermal_postprocessing_path=str(alternative))
    report=dict(thermochemical_gap_validated=True,native_partition_protocol=protocol)
    assert not native_partition_dependencies_ok(report,left,right)
    report['native_partition_dependencies']={k:hashlib.sha256(p.read_bytes()).hexdigest()
                                           for k,p in [('chosen',chosen),('alternative',alternative)]}
    assert native_partition_dependencies_ok(report,left,right)
    alternative.write_text('{"spin":1}')
    assert not native_partition_dependencies_ok(report,left,right)
    assert native_partition_dependencies_ok({'thermochemical_gap_validated':False},{},{})

