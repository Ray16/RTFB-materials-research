import pytest

from redox.core.common import ROOT
from redox.workflows import polymer_campaign as api


def module():
    return api.load_module('first_site_grafting_tests', ROOT / 'scripts/validation/analyze_first_site_grafting.py')


def test_shared_reagent_shift_cancels_from_local_polymer_difference():
    m = module()
    first = m.deltas(-100., -70., -60., -31., {'free_ions': 28., 'contact_pair': 27.9})
    shifted = m.deltas(-100., -70., -60., -31., {'free_ions': 128., 'contact_pair': 127.9})
    assert first['polymer_minus_capped_kJ_mol'] == shifted['polymer_minus_capped_kJ_mol']
    for name in first['first_graft_kJ_mol']:
        assert shifted['first_graft_kJ_mol'][name] - first['first_graft_kJ_mol'][name] == pytest.approx(100 * m.EV_KJ)


@pytest.mark.parametrize('value', [None, float('nan'), float('inf')])
def test_missing_or_invalid_G_cannot_be_replaced_by_electronic_energy(value):
    m = module()
    with pytest.raises((TypeError, ValueError)):
        m.deltas(value, -70., -60., -31., {'free_ions': 28.})


def test_element_and_charge_imbalances_fail_independently():
    m = module()
    with pytest.raises(ValueError, match='atom/charge'):
        m.balance(['CC'], ['CO'])
    with pytest.raises(ValueError, match='atom/charge'):
        m.balance(['[Cl-]'], ['[Cl]'])
