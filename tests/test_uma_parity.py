"""Incomplete states and matching energies with incorrect forces cannot pass."""
from copy import deepcopy

import pytest

from redox.validation.uma_parity import compare_points


def records():
    identity = dict(model='test', checkpoint_sha256='weights', reference_sha256={'atom_refs':'atoms', 'form_elem_refs':'elements'}, task_name='omol')
    case = dict(id='water', charge=0, mult=1, atoms=3, frames=['frame.xyz'])
    plan = dict(model_identity=identity, cases=[case], minimum_energy_tolerance_eV=.002,
        energy_tolerance_eV_per_atom=1e-4, force_RMSE_tolerance_eV_A=.002, force_max_tolerance_eV_A=.01)
    point = dict(status='UMA_point_evaluations_complete_cross_mode_review_pending', model_identity=identity,
        cases=[dict(case, frames=[dict(energy_eV=-2., forces_eV_A=[[0.,0.,0.]]*3)])])
    return plan, point, deepcopy(point)


def test_complete_identical_points_pass_but_do_not_promote_production():
    report = compare_points(*records())
    assert report['all_energy_and_force_checks_pass']
    assert not report['production_migration_qualified']


def test_matching_energies_do_not_mask_bad_forces():
    plan, a, b = records()
    b['cases'][0]['frames'][0]['forces_eV_A'] = [[.1,0.,0.]]*3
    report = compare_points(plan, a, b)
    assert not report['all_energy_and_force_checks_pass']
    assert report['comparisons'][0]['frames'][0]['checks']['energy']


@pytest.mark.parametrize('change', ['missing_case', 'wrong_spin', 'reference', 'nan', 'missing_frame'])
def test_invalid_or_incomplete_comparisons_are_refused(change):
    plan, a, b = records()
    if change == 'missing_case': b['cases'] = []
    elif change == 'wrong_spin': b['cases'][0]['mult'] = 3
    elif change == 'reference': b['model_identity'] = dict(b['model_identity'], reference_sha256={})
    elif change == 'nan': b['cases'][0]['frames'][0]['energy_eV'] = float('nan')
    else: b['cases'][0]['frames'] = []
    with pytest.raises(ValueError): compare_points(plan, a, b)


def test_diagnostic_behavior_override_remains_diagnostic_when_numbers_match():
    plan, a, b = records()
    b['diagnostic_only'] = True
    report = compare_points(plan, a, b)
    assert report['all_energy_and_force_checks_pass'] and report['diagnostic_only']
    assert not report['production_migration_qualified']
