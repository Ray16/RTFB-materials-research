"""Sparse or phase-inconsistent data cannot tighten the reported class prior."""
import json

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
