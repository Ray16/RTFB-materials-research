from itertools import combinations

import pytest
from rdkit import Chem
from rdkit.Chem import Descriptors

from redox.build.polymer_cases import build_case
from redox.properties.polymer_capacity import fragment_capacity_bound, target_exceeds_upper_bound


FRAGMENT = "[*:1]n1c(=O)c2cc3c(=O)n(CC)c(=O)c3cc2c1=O"


def test_analytic_mass_matches_all_eight_exact_occupancy_graphs():
    for n in range(4):
        masses = []
        for occupied in combinations([2, 3, 4], n):
            case = build_case(FRAGMENT, 5, [2, 3, 4], list(occupied), ["CW"] * 4)
            graph_mass = Descriptors.MolWt(Chem.MolFromSmiles(case["smiles"]))
            bound = fragment_capacity_bound(FRAGMENT, 5, 3, n, 2)
            assert bound["mass_g_mol"] == pytest.approx(graph_mass, abs=1e-8)
            masses.append(graph_mass)
        assert max(masses) - min(masses) < 1e-8


def test_more_grafts_increase_the_bound_but_inactive_mass_reduces_it():
    values = [fragment_capacity_bound(FRAGMENT, 5, 3, n, 2)["capacity_upper_bound_mAh_g"] for n in range(4)]
    assert values == sorted(values) and values[0] == 0
    plain = fragment_capacity_bound(FRAGMENT, 5, 3, 2, 2)
    diluted = fragment_capacity_bound(FRAGMENT, 5, 3, 2, 2, additional_inactive_mass_g_mol=100)
    assert diluted["capacity_upper_bound_mAh_g"] < plain["capacity_upper_bound_mAh_g"]


def test_electron_budget_is_an_explicit_assumption_not_accessibility():
    one = fragment_capacity_bound(FRAGMENT, 5, 3, 2, 1)
    two = fragment_capacity_bound(FRAGMENT, 5, 3, 2, 2)
    assert two["capacity_upper_bound_mAh_g"] == pytest.approx(2 * one["capacity_upper_bound_mAh_g"])
    assert not two["accessible_capacity_established"]
    assert two["capacity_lower_bound_mAh_g"] == 0
    assert target_exceeds_upper_bound(two, two["capacity_upper_bound_mAh_g"] + 1)
    assert not target_exceeds_upper_bound(two, 1)


@pytest.mark.parametrize("args", [(5, 3, 4, 2), (5, 6, 1, 2), (True, 1, 1, 2), (5, 3, 1, 0)])
def test_invalid_compositions_or_electron_counts_are_rejected(args):
    with pytest.raises(ValueError):
        fragment_capacity_bound(FRAGMENT, *args)


def test_no_handles_has_a_defined_zero_bound_without_a_conversion_claim():
    result = fragment_capacity_bound(FRAGMENT, 5, 0, 0, 2)
    assert result["capacity_upper_bound_mAh_g"] == 0
    assert result["reactive_handle_conversion"] is None


def test_viologen_route_cannot_silently_drop_its_counterion_mass():
    charged = "[*:1][n+]1ccc(-c2cc[n+](C)cc2)cc1"
    with pytest.raises(ValueError, match="ionic mass convention"):
        fragment_capacity_bound(charged, 5, 3, 2, 2)
