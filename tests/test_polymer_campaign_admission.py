"""Prevent foreign cases and electronic-only records from becoming native deltaG."""

import json
import pytest
from redox.workflows import polymer_campaign as campaign


@pytest.mark.parametrize("problem", ["wrong_site", "electronic_only", "old_thermal"])
def test_native_G_requires_exact_task_identity_and_qualified_thermal_definition(
    tmp_path, monkeypatch, problem
):
    task = {"key": "case", "occupied_sites": [2, 3], "smiles": "CC"}
    (tmp_path / "plan.json").write_text(json.dumps({"tasks": [task]}))
    folder = tmp_path / "thermal/done"
    folder.mkdir(parents=True)
    (folder / "case.json").write_text("{}")
    report = dict(
        task=dict(task),
        G_eV=-100.0,
        E_smd_eV=-101.0,
        thermal_QC="ok",
        thermal_definition="native_solution_dft_qrrho_spin",
    )
    if problem == "wrong_site":
        report["task"]["occupied_sites"] = [2, 4]
    elif problem == "electronic_only":
        del report["G_eV"]
    else:
        report["thermal_definition"] = "legacy_gfn2_gas_rrho"
    monkeypatch.setattr(campaign, "verify_state_checkpoint", lambda path: report)
    records, invalid = campaign.admitted_records(tmp_path, "thermal")
    assert not records and len(invalid) == 1


def test_native_reagent_terms_stay_blank_when_reference_receipts_are_missing(
    tmp_path, monkeypatch
):
    (tmp_path / "cases.json").write_text(json.dumps({"fragment": "[*:1]Oc1ccccc1"}))

    def missing(path):
        raise FileNotFoundError(str(path))

    monkeypatch.setattr(campaign, "verify_state_checkpoint", missing)
    terms = campaign.reaction_terms(tmp_path, "new_phenol", "dbu")
    assert terms["constant_eV"] == {"free_ions": None, "contact_pair": None}
    assert len(terms["missing_reagents"]) == 5


def test_electronic_receipt_without_independent_geometry_review_is_rejected(
    tmp_path, monkeypatch
):
    task = {"key": "case", "smiles": "CC"}
    (tmp_path / "plan.json").write_text(json.dumps({"tasks": [task]}))
    folder = tmp_path / "energy/done"
    folder.mkdir(parents=True)
    (folder / "case.json").write_text("{}")
    monkeypatch.setattr(
        campaign,
        "verify_state_checkpoint",
        lambda path: {"task": task, "E_smd_eV": -100.0},
    )
    records, invalid = campaign.admitted_records(tmp_path, "energy")
    assert not records and len(invalid) == 1
