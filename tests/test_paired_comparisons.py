"""An incomplete control must not manufacture a matched substituent effect."""
import importlib.util
from pathlib import Path

import pytest

from redox.screening.comparisons import paired_potential_shifts


def row(gid, **kwargs):
    return dict(id=gid, registry_id=gid,
                states='[["rest",2,1,0],["first",1,2,-1],["second",0,1,-2]]',
                **kwargs)


def test_partial_control_does_not_renumber_or_use_surface_estimate():
    control = row('base', E1_vs_Fc_V=None, E2_vs_Fc_V=-1.0, E1_surface_estimate_V=-.6)
    candidate = row('new', comparison_baseline='base', E1_vs_Fc_V=-.3, E2_vs_Fc_V=-.8)
    result, = paired_potential_shifts([control, candidate])
    assert result['delta_E1_V'] is None
    assert result['delta_E2_V'] == pytest.approx(.2)
    assert result['status'] == 'one potential shift admitted'


def test_signed_shifts_and_reference_free_spacing_are_retained():
    control = row('base', E1_vs_Fc_V=-.5, E2_vs_Fc_V=-1., wave_spacing_V=.5)
    candidate = row('new', comparison_baseline='base', E1_vs_Fc_V=-.7,
                    E2_vs_Fc_V=-.9, wave_spacing_V=.2)
    result, = paired_potential_shifts([control, candidate])
    assert result['delta_E1_V'] == pytest.approx(-.2)
    assert result['delta_E2_V'] == pytest.approx(.1)
    assert result['delta_wave_spacing_V'] == pytest.approx(-.3)
    assert not result['experimental_improvement_verified']


def test_missing_control_and_incompatible_charge_paths():
    candidate = row('new', comparison_baseline='base', E1_vs_Fc_V=-.3)
    assert paired_potential_shifts([candidate])[0]['status'] == 'pending control'
    with pytest.raises(ValueError, match='Incompatible'):
        paired_potential_shifts([candidate, row('base') | {'states': [['neu', 0, 1, 0], ['red', -1, 2, -1]]}])


def test_older_configuration_control_is_collected_once():
    spec = importlib.util.spec_from_file_location('control_collection',
        Path(__file__).resolve().parents[1] / 'scripts/pipeline/collect_verified_candidates.py')
    module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    base = dict(id='older_control')
    candidates = [dict(id='new', registry_id='new', comparison_baseline='older_control')]
    assert module.comparison_control_groups(candidates, {'key': base}) == [base]
    assert module.comparison_control_groups(candidates + [dict(registry_id='older_control')], {'key': base}) == []
    with pytest.raises(ValueError, match='Unregistered comparison controls'):
        module.comparison_control_groups(candidates, {})
