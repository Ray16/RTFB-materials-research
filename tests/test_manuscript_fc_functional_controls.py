"""Reject reference corrections or numerical claims unsupported by the raw states."""
from copy import deepcopy
import importlib.util
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).parents[1]
spec = importlib.util.spec_from_file_location('fc_claims_test', ROOT / 'manuscript/check_claims.py')
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


@pytest.fixture
def data():
    folder = ROOT / 'manuscript/evidence/2026-10-25_reference_diagnostic'
    names = ['Fc_functional_review.json', 'benchmark.json'] + [
        f'Fc_functional_{method}_{state}.json'
        for method in ('wb97mv', 'b3lyp', 'pbe0') for state in ('neu', 'ox')]
    return {name: json.loads((folder / name).read_text()) for name in names}


def test_six_electronic_states_reconstruct_reference_diagnostics(data):
    result = module.audit_fc_functional_controls(data.__getitem__)
    assert result['electronic_states'] == 6
    assert result['gaps_and_counterfactuals_reconstructed']
    assert not result['alternative_potential_qualified']


@pytest.mark.parametrize('change', ['energy', 'counterfactual', 'spin', 'nan', 'state', 'promotion'])
def test_unsupported_Fc_claims_are_rejected(data, change):
    data = deepcopy(data)
    row = data['Fc_functional_b3lyp_ox.json']
    if change == 'energy':
        row['e_smd_eV'] += .1
    elif change == 'counterfactual':
        data['Fc_functional_review.json']['methods']['b3lyp']['reference_only_counterfactual']['MAE_V'] = .1
    elif change == 'spin':
        row['smd_spin_contam'] = .2
    elif change == 'nan':
        row['e_smd_eV'] = float('nan')
    elif change == 'state':
        row['charge'] = 0
    else:
        data['Fc_functional_review.json']['alternative_potentials_calculated'] = True
    with pytest.raises(ValueError):
        module.audit_fc_functional_controls(data.__getitem__)
