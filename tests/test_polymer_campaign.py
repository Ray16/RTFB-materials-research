"""Condition-specific grafting and neighbor penalties cannot be conflated."""

import pytest
from redox.core.common import EV_KJ
from redox.workflows.polymer_campaign import nonadditivity, conditional_grafting


def test_precursor_and_base_terms_cancel_from_neighbor_penalty():
    values = (-100.0, -110.0, -111.0, -119.0)
    a = conditional_grafting(*values, reagent_constant=5.0)
    b = conditional_grafting(*values, reagent_constant=8.0)
    assert a["neighbor_penalty_deltaDeltaG_kJ_mol"] == pytest.approx(2 * EV_KJ)
    assert b["neighbor_penalty_deltaDeltaG_kJ_mol"] == pytest.approx(
        a["neighbor_penalty_deltaDeltaG_kJ_mol"]
    )
    assert b["graft_j_given_i_kJ_mol"] - a["graft_j_given_i_kJ_mol"] == pytest.approx(
        3 * EV_KJ
    )
    assert a["graft_j_on_empty_kJ_mol"] < a["graft_j_given_i_kJ_mol"] < 0


def test_reverse_order_has_the_same_nonadditivity_with_site_specific_first_grafting():
    result = conditional_grafting(-100.0, -110.0, -111.0, -119.0, 5.0)
    assert result["graft_j_given_i_kJ_mol"] - result[
        "graft_j_on_empty_kJ_mol"
    ] == pytest.approx(result["neighbor_penalty_deltaDeltaG_kJ_mol"])
    assert result["graft_i_given_j_kJ_mol"] - result[
        "graft_i_on_empty_kJ_mol"
    ] == pytest.approx(result["neighbor_penalty_deltaDeltaG_kJ_mol"])
    assert result["graft_i_on_empty_kJ_mol"] != result["graft_j_on_empty_kJ_mol"]


@pytest.mark.parametrize("bad", [None, float("nan"), float("inf"), True])
def test_missing_or_invalid_values_are_not_promoted_to_grafting_results(bad):
    with pytest.raises(ValueError):
        nonadditivity(-100.0, bad, -111.0, -119.0)
    with pytest.raises(ValueError):
        conditional_grafting(-100.0, -110.0, -111.0, -119.0, bad)
