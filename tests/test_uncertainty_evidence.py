"""Sparse or phase-inconsistent data cannot tighten the reported class prior."""
import json
import csv

import pytest

from redox.core import common
from redox.core.protocol import ACTIVE_SP
from redox.properties import potentials
from redox.screening import scorecard


def test_charge_class_fallback_requires_multiple_current_benchmark_events(tmp_path, monkeypatch):
    cfg = common.load_config('scorecard_config')
    row = dict(status='ok', calibration_eligible='True', energy_protocol=json.dumps(ACTIVE_SP),
               reference_abs_V=4., reference_thermal_definition='native_solution_dft_qrrho_spin',
               thermal_definition='native_solution_dft_qrrho_spin', charge_ox=2, calc=-1., err=.1)
    benchmark = [row]
    monkeypatch.setattr(scorecard, 'RESULTS', tmp_path)
    monkeypatch.setattr(scorecard, '_load', lambda name: benchmark if name=='orop_benchmark' else [])
    monkeypatch.setattr(scorecard, 'experimental_e_points', lambda rows: [dict(family='q',res=.05)])
    monkeypatch.setattr(potentials, 'fc_live', lambda: 4.)
    monkeypatch.setattr(common, 'read_result', lambda *a: dict(thermal_definition='native_solution_dft_qrrho_spin'))
    uncertainty = scorecard.derive_uncertainties([])
    assert uncertainty['pooled_redox'] is None
    assert uncertainty['redox_calibration_status'] == 'prior_limited'
    assert uncertainty['charge_redox']['multi'] == cfg.SIGMA_REDOX_V['multi']
    assert cfg.sigma_redox(2,1,charge_sigma=uncertainty['charge_redox']) == .8
    benchmark.append(dict(row, err=.2))
    uncertainty = scorecard.derive_uncertainties([])
    assert uncertainty['charge_redox']['multi'] == pytest.approx(.158)
    assert uncertainty['charge_redox']['anion'] == cfg.SIGMA_REDOX_V['anion']


@pytest.mark.parametrize('pairs,expected,status', [
    ([(.5,.1)], .1, 'prior_limited'),
    ([(.5,.1),(.7,.3)], .1, 'computational_reference_only'),
    ([(.2,.1),(.8,.1),(1.2,.1)], .593, 'computational_reference_only'),
])
def test_computational_lambda_spread_cannot_imply_empirical_or_sparse_precision(
        tmp_path, monkeypatch, pairs, expected, status):
    directory = tmp_path / 'd3tales_reorg_validation_prod'
    directory.mkdir()
    with (directory / 'comparison.csv').open('w', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(['our_smd_eV', 'd3tales_eV'])
        writer.writerows(pairs)
    monkeypatch.setattr(scorecard, 'RESULTS', tmp_path)
    monkeypatch.setattr(scorecard, '_load', lambda name: [])
    monkeypatch.setattr(scorecard, 'experimental_e_points', lambda rows: [])
    monkeypatch.setattr(potentials, 'fc_live', lambda: 4.)
    monkeypatch.setattr(common, 'read_result', lambda *a: dict(thermal_definition='native_solution_dft_qrrho_spin'))
    value = scorecard.derive_uncertainties([])
    assert value['sigma_lambda'] == pytest.approx(expected)
    assert value['lambda_comparison_events'] == len(pairs)
    assert value['lambda_calibration_status'] == status
    assert value['lambda_empirical_solution_accuracy_validated'] is False
    assert value['lambda_uncertainty_kind'] == 'computational_model_discrepancy_scale_with_prior_floor'
    if len(pairs) == 2:
        # Identical model discrepancies have zero spread, but a nonzero bias.
        assert value['lambda_comparison_spread_eV'] == pytest.approx(0.)
        assert value['lambda_comparison_median_residual_eV'] == pytest.approx(.4)
