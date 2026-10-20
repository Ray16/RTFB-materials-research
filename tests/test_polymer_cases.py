"""Independent topology, stereo and stoichiometry checks for occupancy scenarios."""

import json
from pathlib import Path
import pytest
from rdkit import Chem
from redox.build.polymer_cases import (
    build_case,
    panel,
    stereo_fingerprint,
    resolve_motif,
)

SPEC = json.loads(
    (
        Path(__file__).resolve().parents[1] / "config/polymer_spacing_pilot.json"
    ).read_text()
)
FRAGMENT = "[*:1]N1C(=O)c2ccccc2C1=O"


def test_matched_panel_keeps_backbone_stereo_and_distinct_site_chemistry():
    cases, cycles = panel(SPEC, FRAGMENT)
    assert len(cases) == 6 and len(cycles) == 2
    for case in cases:
        assert case["backbone_carbon_count"] == 10 and case["backbone_CC_bonds"] == 9
        assert case["site_states"][0] == case["site_states"][-1] == "styrene"
        assert case["site_states"].count("grafted") == len(case["occupied_sites"])
        mol = Chem.MolFromSmiles(case["smiles"])
        assert sum(a.GetAtomicNum() == 17 for a in mol.GetAtoms()) == 3 - len(
            case["occupied_sites"]
        )
        assert len(Chem.GetMolFrags(mol)) == 1
        # Canonicalization may reorder neighbors but must preserve orientation.
        remapped = Chem.MolToSmiles(
            Chem.MolFromSmiles(case["mapped_smiles"]),
            canonical=False,
            isomericSmiles=True,
        )
        assert stereo_fingerprint(remapped, 5) == case["backbone_stereo_fingerprint"]
    assert len({str(c["backbone_stereo_fingerprint"]) for c in cases}) == 1


def test_adjacent_and_spaced_pairs_are_composition_matched():
    cases, cycles = panel(SPEC, FRAGMENT)
    lookup = {c["key"]: c for c in cases}
    a, b = lookup["occ_2_3"], lookup["occ_2_4"]
    assert (
        a["atoms"] == b["atoms"]
        and a["model_loading_per_repeat"] == b["model_loading_per_repeat"] == 0.4
    )
    assert [c["intervening_repeat_units"] for c in cycles] == [0, 1]


def test_unreactive_styrene_is_not_silently_changed_to_benzyl_chloride():
    plain = build_case(FRAGMENT, 3, [], [], ["CW", "CW"])
    assert plain["site_states"] == ["styrene"] * 3
    assert not any(
        a.GetAtomicNum() == 17 for a in Chem.MolFromSmiles(plain["smiles"]).GetAtoms()
    )


def test_grafting_a_site_without_a_handle_is_rejected():
    with pytest.raises(ValueError, match="reactive"):
        build_case(FRAGMENT, 5, [2, 3, 4], [1, 2], ["CW"] * 4)


@pytest.mark.parametrize(
    "motif_id",
    [
        "exp_phthalimide_control",
        "design_pmdi_ethyl",
        "design_ndi_ethyl",
        "exp_bq_oxy_methoxy",
        "exp_aq_1oxy",
    ],
)
def test_the_same_panel_applies_to_N_and_O_grafted_motif_families(motif_id):
    root = Path(__file__).resolve().parents[1]
    feed = json.loads((root / "results/candidates/current.json").read_text())
    row = next(
        row
        for row in feed["original39"] + feed["additional_candidates"]
        if row["id"] == motif_id
    )
    cases, cycles = panel(SPEC, row["frag"])
    assert len(cases) == 6 and len(cycles) == 2
    assert all(case["charge"] == 0 and case["chain_ends"] == 2 for case in cases)


def test_charged_motifs_require_a_separate_ionic_model_before_submission():
    with pytest.raises(ValueError, match="counterion"):
        panel(SPEC, "[*:1][n+]1ccccc1")


def test_new_neutral_motif_can_be_defined_without_candidate_registry():
    spec = {
        "motif": {
            "id": "custom_phenol",
            "fragment": "[*:1]Oc1ccccc1",
            "precursor_smiles": "Oc1ccccc1",
        }
    }
    motif, source = resolve_motif(spec, {})
    assert motif["id"] == "custom_phenol" and source == "inline configuration"
    cases, cycles = panel(SPEC, motif["frag"])
    assert len(cases) == 6 and len(cycles) == 2


@pytest.mark.parametrize("smiles", ["CCO", "not_smiles"])
def test_inline_precursor_must_match_the_exact_attachment_fragment(smiles):
    with pytest.raises(ValueError, match="precursor"):
        resolve_motif(
            {
                "motif": {
                    "id": "custom",
                    "fragment": "[*:1]Oc1ccccc1",
                    "precursor_smiles": smiles,
                }
            },
            {},
        )
