"""Join systematic designs to the existing evaluator without inheriting properties.

Exact graft graphs identify candidates; declared controls identify comparisons.
Pool coverage, admitted properties, route admission and polymer evidence remain
separate dimensions. This consumer never registers a candidate or launches QM.
"""

from collections import Counter
import math

from rdkit import Chem

from redox.screening.comparisons import paired_potential_shifts


def identity(smiles):
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise ValueError("Invalid exact graft SMILES")
    if len(Chem.GetMolFrags(molecule)) != 1:
        raise ValueError("Specify a single exact graft; ionic models stay separate")
    return Chem.MolToSmiles(molecule, isomericSmiles=True), Chem.MolToInchiKey(molecule)


def finite(value):
    return type(value) in (int, float) and math.isfinite(value)


def attach_evaluations(pool, feed):
    """Canonical graph join; retain aliases and all symmetry-equivalent recipes."""
    records = {}
    ids = {}
    for item in pool:
        smiles, key = identity(item["grafted_smiles"])
        if key != item["inchikey"]:
            raise ValueError("Pool identity does not match its declared graph")
        if key in records and records[key]["grafted_smiles"] != smiles:
            raise ValueError("InChIKey collision between distinct exact graphs")
        record = records.setdefault(
            key,
            dict(
                inchikey=key,
                grafted_smiles=smiles,
                core=item["core"],
                pool_ids=[],
                recipes=[],
                evaluation=None,
                evaluation_role=None,
                fragment=item["group"]["frag"],
                states=item["group"]["states"],
            ),
        )
        record["pool_ids"].append(item["group"]["id"])
        for recipe in item["recipes"]:
            if recipe not in record["recipes"]:
                record["recipes"].append(recipe)
    for role, rows in (
        ("candidate", feed["original39"] + feed["additional_candidates"]),
        ("comparison_control", feed.get("comparison_controls", [])),
    ):
        for row in rows:
            smiles, key = identity(row["exact_grafted_smiles"])
            if key != row["exact_grafted_inchikey"]:
                raise ValueError("Evaluation identity does not match its graph")
            identifier = row.get("registry_id") or row["id"]
            if identifier in ids and ids[identifier] != key:
                raise ValueError("Candidate ID belongs to different exact graphs")
            ids[identifier] = key
            record = records.setdefault(
                key,
                dict(
                    inchikey=key,
                    grafted_smiles=smiles,
                    core=row.get("search_core") or row["family"],
                    pool_ids=[],
                    recipes=[],
                    evaluation=None,
                    evaluation_role=None,
                    fragment=row["frag"],
                    states=row["states"],
                ),
            )
            if record["grafted_smiles"] != smiles:
                raise ValueError("Evaluation graph differs from the pool identity")
            if record["evaluation"] is not None:
                previous = record["evaluation"]
                if previous != row:
                    raise ValueError("Ambiguous duplicate exact-graph evaluation")
                continue
            if row.get("parent_metrics_inherited") is True:
                raise ValueError(
                    "Parent measurements cannot enter exact-graft evaluation"
                )
            for field in (
                "E1_vs_Fc_V",
                "E2_vs_Fc_V",
                "wave_spacing_V",
                "MW",
                "accessible_capacity_mAh_g",
            ):
                if row.get(field) is not None and not finite(row[field]):
                    raise ValueError(
                        f"Nonfinite or nonnumeric evaluator field: {field}"
                    )
            record.update(evaluation=row, evaluation_role=role, fragment=row["frag"])
    return sorted(records.values(), key=lambda record: record["inchikey"])


def exact_evidence(row, evidence, product_field):
    """Never attach energy evidence merely because an ID happens to match."""
    if evidence is not None:
        product = evidence.get(product_field)
        if (
            product is None
            or identity(product)[0] != identity(row["exact_grafted_smiles"])[0]
        ):
            raise ValueError("Grafting evidence belongs to a different exact product")
    return evidence


def build_review(pool, feed, grafting, bases, assessment, polymers, plan, screen):
    records = attach_evaluations(pool, feed)
    evaluated = (
        feed["original39"]
        + feed["additional_candidates"]
        + feed.get("comparison_controls", [])
    )
    comparisons = paired_potential_shifts(evaluated)
    by_id = {row.get("registry_id") or row["id"]: row for row in evaluated}
    graft_by_id = {row["id"]: row for row in grafting["candidates"]}
    base_by_id = {row["id"]: row for row in bases["candidates"]}
    polymer_by_fragment = {}
    for model in polymers["models"]:
        key = Chem.MolToSmiles(Chem.MolFromSmiles(model["fragment"]))
        polymer_by_fragment.setdefault(key, []).append(model)
    for record in records:
        row = record.pop("evaluation")
        identifier = (
            (row.get("registry_id") or row["id"]) if row else record["pool_ids"][0]
        )
        group = row or dict(id=identifier, frag=record["fragment"])
        route = screen(group)
        pair = bool(
            row and finite(row.get("E1_vs_Fc_V")) and finite(row.get("E2_vs_Fc_V"))
        )
        selected = assessment.get("candidates", {}).get(identifier, {})
        reference = (
            exact_evidence(row, graft_by_id.get(identifier), "product_smiles")
            if row
            else None
        )
        base = base_by_id.get(identifier) if row else None
        if base:
            exact_evidence(row, base["reaction"], "product_smiles")
        fragment = Chem.MolToSmiles(Chem.MolFromSmiles(record["fragment"]))
        model_evidence = polymer_by_fragment.get(fragment, [])
        checked_models = [
            model
            for model in model_evidence
            if model.get("stage") == "stationary_DFT_electronic_structure"
            and model.get("structural_review", {}).get("issues") == []
        ]
        record.update(
            id=identifier,
            name=row["name"] if row else identifier,
            precursor_smiles=route.get("precursor_smiles"),
            reaction_class=route.get("reaction_class"),
            attachment_definition=(
                "The single mapped dummy in fragment defines the benzyl-carbon bond"
            ),
            source_urls=row.get("source_urls") if row else None,
            design_origin="systematic_enumeration"
            if record["pool_ids"]
            else "source_guided_or_earlier_cohort",
            evaluation_stage="admitted_potential_pair"
            if pair
            else "partial_or_QC_pending"
            if row
            else "proposed_only",
            E1_vs_Fc_V=row.get("E1_vs_Fc_V") if row else None,
            E2_vs_Fc_V=row.get("E2_vs_Fc_V") if row else None,
            wave_spacing_V=row.get("wave_spacing_V") if row else None,
            MW=row.get("MW") if row else None,
            accessible_capacity_mAh_g=row.get("accessible_capacity_mAh_g")
            if row
            else None,
            redox_QC=row.get("solution_potentials_status") if row else "not_computed",
            model_validation_status=row.get("model_validation_status") if row else None,
            remaining_work=row.get("remaining_work")
            if row
            else "Route review before full redox computation",
            baseline=row.get("comparison_baseline") if row else None,
            design_axis=(
                row.get("search_contrast") or plan["contrast_overrides"].get(identifier)
            )
            if row
            else None,
            route_status=route["route_status"],
            prospective_full_redox_admitted=route["proceed_to_full_redox"],
            missing_route_evidence=route.get("missing_route_evidence", []),
            selection_label=selected.get("display_label"),
            reference_grafting_deltaG_kJ_mol=reference.get("deltaG_reference_kJ_mol")
            if reference
            else None,
            reference_grafting_status=reference.get("status")
            if reference
            else "not_evaluated",
            base_models=base.get("base_models") if base else None,
            checked_polymer_models=len(checked_models),
            polymer_model_keys=[model["key"] for model in checked_models],
            polymer_scope=(
                "Exact attachment-fragment method controls; not capped-site redox "
                "or whole-polymer properties"
            ),
            synthesis_eligible=selected.get("eligible_for_synthesis_shortlist", False),
        )
        record["design_question"] = plan["axes"].get(
            record["design_axis"],
            "Enumerated cell; no prospective batch selection declared"
            if not row
            else "Source-guided or earlier cohort; no single-axis hypothesis declared",
        )
        if not row:
            action = "review_route_before_full_redox"
        elif row.get("model_primary_admissible") is False:
            action = "resolve_scoped_model_validation_limits"
        elif not pair:
            action = "resolve_remaining_redox_fields_and_QC"
        elif record["baseline"] and not all(
            finite(by_id.get(record["baseline"], {}).get(k))
            for k in ("E1_vs_Fc_V", "E2_vs_Fc_V")
        ):
            action = "complete_declared_control_before_comparing"
        else:
            action = "review_matched_effect_and_polymer_environment"
        record["evaluation_feedback"] = action
        # Route admission stays independent of historical property completion.
        record["prospective_route_action"] = (
            "eligible_for_reviewed_campaign"
            if route["proceed_to_full_redox"]
            else "exclude_this_route"
            if route["route_status"].startswith("exclude")
            else "resolve_candidate_specific_route_evidence"
        )
    enriched = []
    for pair in comparisons:
        row = by_id[pair["id"]]
        baseline = by_id.get(pair["baseline"])
        axis = row.get("search_contrast") or plan["contrast_overrides"].get(pair["id"])
        enriched.append(
            dict(
                pair,
                core=row.get("search_core") or row["family"],
                contrast=axis,
                design_question=plan["axes"].get(
                    axis, "Declared source-guided comparison; inspect scope"
                ),
                delta_MW=row["MW"] - baseline["MW"]
                if baseline and finite(row.get("MW")) and finite(baseline.get("MW"))
                else None,
                next_question=(
                    "Inspect both wave shifts, mass cost, route and polymer behavior; "
                    "no automatic improvement claim"
                ),
            )
        )
    coverage = []
    for core in sorted({record["core"] for record in records}):
        members = [record for record in records if record["core"] == core]
        coverage.append(
            dict(
                family=core,
                unique_designs=len(members),
                in_enumerated_pool=sum(bool(row["pool_ids"]) for row in members),
                selected_candidates=sum(
                    row["evaluation_role"] == "candidate" for row in members
                ),
                external_controls=sum(
                    row["evaluation_role"] == "comparison_control" for row in members
                ),
                admitted_potential_pairs=sum(
                    row["evaluation_stage"] == "admitted_potential_pair"
                    for row in members
                ),
                prospective_route_proceed=sum(
                    row["prospective_full_redox_admitted"] for row in members
                ),
            )
        )
    cells = {}
    for record in records:
        for recipe in record["recipes"]:
            key = tuple(recipe[k] for k in ("core", "anchor", "cap", "ring"))
            cells.setdefault(key, set()).add(record["inchikey"])
    design_cells = [
        dict(
            zip(("core", "anchor", "cap", "ring"), key),
            unique_graphs=len(values),
            selected_graphs=sum(
                row["inchikey"] in values and row["evaluation_role"] == "candidate"
                for row in records
            ),
        )
        for key, values in sorted(cells.items())
    ]
    return dict(
        schema=1,
        strategy=plan["strategy"],
        source_candidate_as_of=feed.get("as_of_Chicago"),
        records=records,
        comparisons=enriched,
        family_coverage=coverage,
        design_cells=design_cells,
        counts=dict(
            unique_designs=len(records),
            enumerated_pool=len({identity(row["grafted_smiles"])[1] for row in pool}),
            selected_candidates=len(feed["original39"])
            + len(feed["additional_candidates"]),
            external_controls=len(feed.get("comparison_controls", [])),
            evaluation_stages=dict(Counter(row["evaluation_stage"] for row in records)),
            prospective_route_proceed=sum(
                row["prospective_full_redox_admitted"] for row in records
            ),
            matched_comparisons=len(comparisons),
            complete_matched_pairs=sum(
                row["status"] == "both potential shifts admitted" for row in comparisons
            ),
        ),
        future_axes=plan["future_axes"],
        prospective_rules=plan["prospective_rules"],
        limits=[
            "Enumeration does not establish synthesis feasibility or usable electrons",
            "Historical results remain separate from prospective route admission",
            "Parent measurements and surface estimates never fill missing graft fields",
            "Segment controls do not establish motif-only binding, coupling or rates",
            "No composite score, automatic shortlist or unreviewed new submission",
        ],
    )
