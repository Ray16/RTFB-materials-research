"""The design/evaluation join must preserve identity, missing fields and route scope."""

from copy import deepcopy

import pytest

from redox.screening.expansion_review import attach_evaluations, build_review, identity


def candidate(identifier="a", smiles="CC", **values):
    canonical, key = identity(smiles)
    row = dict(
        id=identifier,
        registry_id=identifier,
        name=identifier,
        exact_grafted_smiles=canonical,
        exact_grafted_inchikey=key,
        family="test",
        frag="[*:1]O",
        states=[["neu", 0, 1, 0], ["red", -1, 2, -1]],
        E1_vs_Fc_V=None,
        E2_vs_Fc_V=None,
    )
    row.update(values)
    return row


def pool_item(identifier="proposed_a", smiles="CC"):
    canonical, key = identity(smiles)
    return dict(
        grafted_smiles=canonical,
        inchikey=key,
        core="test",
        group=dict(id=identifier, frag="[*:1]O", states=[["neu", 0, 1, 0]]),
        recipes=[
            dict(
                core="test", anchor="direct", cap="ethyl", ring="H", template_atom=None
            )
        ],
    )


def feed(*rows):
    return dict(original39=[], additional_candidates=list(rows), comparison_controls=[])


def screen(row):
    return dict(
        route_status="hold",
        proceed_to_full_redox=False,
        missing_route_evidence=["exact_precursor_access_reviewed"],
    )


def review(pool, data, grafting=None, bases=None, assessment=None):
    plan = dict(
        strategy="matched",
        axes={},
        contrast_overrides={},
        future_axes=[],
        prospective_rules=[],
    )
    return build_review(
        pool,
        data,
        grafting or {"candidates": []},
        bases or {"candidates": []},
        assessment or {},
        {"models": []},
        plan,
        screen,
    )


def test_graph_aliases_join_without_duplicate_candidate_or_lost_recipes():
    first = pool_item()
    second = deepcopy(first)
    second["group"]["id"] = "symmetry_recipe"
    second["recipes"][0]["template_atom"] = 5
    rows = attach_evaluations([first, second], feed(candidate()))
    assert len(rows) == 1 and len(rows[0]["recipes"]) == 2
    assert rows[0]["evaluation"]["id"] == "a"


def test_wrong_pool_or_evaluator_identity_is_refused():
    wrong = pool_item()
    wrong["inchikey"] = identity("CCC")[1]
    with pytest.raises(ValueError, match="identity"):
        attach_evaluations([wrong], feed())
    wrong = candidate()
    wrong["exact_grafted_inchikey"] = identity("CCC")[1]
    with pytest.raises(ValueError, match="identity"):
        attach_evaluations([], feed(wrong))


def test_same_id_for_different_products_cannot_transfer_results():
    with pytest.raises(ValueError, match="Candidate ID"):
        attach_evaluations([], feed(candidate(), candidate(smiles="CCC")))


def test_ambiguous_duplicate_graph_evaluations_are_not_arbitrarily_chosen():
    with pytest.raises(ValueError, match="Ambiguous"):
        attach_evaluations([], feed(candidate(), candidate(identifier="alias")))


def test_parent_metrics_cannot_enter_exact_graft_ledger():
    with pytest.raises(ValueError, match="Parent"):
        attach_evaluations([], feed(candidate(parent_metrics_inherited=True)))


@pytest.mark.parametrize("value", [float("nan"), float("inf"), True])
def test_published_numeric_fields_cannot_contain_nonfinite_or_boolean_values(value):
    with pytest.raises(ValueError, match="evaluator field"):
        attach_evaluations([], feed(candidate(E1_vs_Fc_V=value)))


def test_proposed_graph_does_not_inherit_evaluation_or_get_synthesis_admission():
    data = review([pool_item(smiles="CCC")], feed(candidate()))
    proposed = next(
        row for row in data["records"] if row["evaluation_stage"] == "proposed_only"
    )
    assert proposed["E1_vs_Fc_V"] is None and proposed["base_models"] is None
    assert (
        not proposed["synthesis_eligible"]
        and not proposed["prospective_full_redox_admitted"]
    )


def test_graft_energy_requires_matching_product_even_when_id_matches():
    with pytest.raises(ValueError, match="different exact product"):
        review(
            [],
            feed(candidate()),
            grafting={"candidates": [dict(id="a", product_smiles="CCC")]},
        )


def test_base_required_candidate_is_retained_without_excluding_positive_HCl_reference():
    result = review(
        [],
        feed(candidate()),
        grafting={
            "candidates": [
                dict(id="a", product_smiles="CC", deltaG_reference_kJ_mol=50)
            ]
        },
        bases={
            "candidates": [
                dict(
                    id="a",
                    reaction={"product_smiles": "CC"},
                    base_models={"dbu": {"free_ions_kJ_mol": -20}},
                )
            ]
        },
        assessment={
            "candidates": {
                "a": {"display_label": "Base required; retained for evaluation"}
            }
        },
    )
    row = result["records"][0]
    assert row["selection_label"].startswith("Base required")
    assert row["base_models"]["dbu"]["free_ions_kJ_mol"] == -20
    assert row["route_status"] == "hold" and not row["synthesis_eligible"]


def test_external_control_is_counted_separately_from_selected_candidates():
    data = feed(candidate())
    data["comparison_controls"] = [candidate(identifier="control", smiles="CCC")]
    result = review([], data)
    assert result["counts"]["selected_candidates"] == 1
    assert result["counts"]["external_controls"] == 1
    assert sum(row["selected_candidates"] for row in result["family_coverage"]) == 1


def test_matched_effect_uses_declared_control_without_surface_estimate_substitution():
    baseline = candidate(identifier="control", E1_vs_Fc_V=-0.5, E2_vs_Fc_V=-1.0)
    variant = candidate(
        identifier="variant",
        smiles="CCC",
        comparison_baseline="control",
        E1_vs_Fc_V=-0.3,
        E2_surface_estimate_V=-0.8,
    )
    (pair,) = review([], feed(baseline, variant))["comparisons"]
    assert pair["delta_E1_V"] == pytest.approx(0.2)
    assert pair["delta_E2_V"] is None


def test_scoped_model_limits_return_to_validation_before_further_design_claims():
    row = candidate(E1_vs_Fc_V=-0.3, E2_vs_Fc_V=-0.8, model_primary_admissible=False)
    result = review([], feed(row))["records"][0]
    assert result["evaluation_feedback"] == "resolve_scoped_model_validation_limits"
    assert not result["synthesis_eligible"]
