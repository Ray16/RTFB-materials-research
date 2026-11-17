"""Systematic enumeration preserves chemistry, coverage and uncertainty scope."""
import pytest
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

from redox.screening import family_search as S
from redox.screening.discovery import families, identify, rediscover


@pytest.fixture(scope="module")
def space():
    return S.enumerate_space()


def test_space_has_all_declared_families_and_no_duplicate_exact_graphs(space):
    assert set(r["core"] for r in space) == {"viologen", "pmdi", "ndi"}
    assert len({r["inchikey"] for r in space}) == len(space)
    assert sum(len(r["recipes"]) for r in space) > len(space)
    for r in space:
        mol=Chem.MolFromSmiles(r["grafted_smiles"])
        assert len(Chem.GetMolFrags(mol)) == 1
        assert ("viologen" if r["core"] == "viologen" else "bis-imide") in families(mol)
        assert Chem.GetFormalCharge(mol) == r["group"]["states"][0][1]
        assert not any(k in r["group"] for k in ("E1", "solubility", "experimental_potential"))


def test_branching_controls_preserve_formula_for_each_core(space):
    for core in ("viologen", "pmdi", "ndi"):
        controls={r["recipes"][0]["cap"]:r for r in space if r["core"] == core
                  and r["recipes"][0]["anchor"] == "direct_N" and r["recipes"][0]["ring"] == "H"}
        formula=lambda r: rdMolDescriptors.CalcMolFormula(Chem.MolFromSmiles(r["grafted_smiles"]))
        assert formula(controls["npropyl"]) == formula(controls["isopropyl"])
        assert formula(controls["nbutyl"]) == formula(controls["isobutyl"])
        assert controls["npropyl"]["inchikey"] != controls["isopropyl"]["inchikey"]


def test_every_anchor_can_be_regenerated_by_declared_grafting_rules(space):
    controls=[r["group"] for r in space if r["recipes"][0]["cap"] == "ethyl"
              and r["recipes"][0]["ring"] == "H"]
    assert len(controls) == 9
    assert all(r["rediscovered"] for r in rediscover(controls,96485.33212))


def test_enumeration_is_deterministic_for_the_same_declared_space(space):
    assert [(r["inchikey"],r["recipes"]) for r in space] == [
        (r["inchikey"],r["recipes"]) for r in S.enumerate_space()]


def test_charged_viologen_precursor_is_a_valid_search_input():
    rows=identify([dict(source_id="monoquat",smiles="C[n+]1ccc(-c2ccncc2)cc1")],96485.33212)
    assert any(r["status"] == "ok" and r["families"] == "viologen" for r in rows)


def test_counterions_are_not_silently_discarded_during_candidate_generation():
    rows=identify([dict(source_id="salt",smiles="C[n+]1ccc(-c2cc[n+](C)cc2)cc1.[Cl-].[Cl-]")],96485.33212)
    assert len(rows) == 1 and rows[0]["status"].startswith("multiple fragments")
    assert rows[0]["grafted_smiles"] is None


def test_comparison_batch_is_balanced_with_matched_ring_sites_and_controls(space):
    selected=S.select_comparisons(space)
    assert len(selected)==12 and len({s['record']['inchikey'] for s in selected})==12
    for core in ('viologen','pmdi','ndi'):
        members=[s for s in selected if s['record']['core']==core]
        assert {s['contrast'] for s in members}=={'cap_length','ring_F','ring_CN','linker_C2'}
        assert len({s['control']['inchikey'] for s in members})==1
        rings=[s for s in members if s['contrast'].startswith('ring_')]
        assert len({s['recipe']['template_atom'] for s in rings})==1
        assert all(s['record']['inchikey']!=s['control']['inchikey'] for s in members)


def test_a_registered_contrast_is_not_silently_replaced(space):
    first=S.select_comparisons(space)[0]['record']['inchikey']
    with pytest.raises(ValueError,match='No new candidate'):
        S.select_comparisons(space,{first})


def test_registered_ring_site_cannot_fall_through_to_another_isomer(space):
    ring=next(s for s in S.select_comparisons(space) if s['contrast']=='ring_F')
    with pytest.raises(ValueError,match='declared site already registered'):
        S.select_comparisons(space,{ring['record']['inchikey']})


def test_donor_panel_keeps_initial_sites_and_controls_after_initial_registration(space):
    initial=S.select_comparisons(space)
    donors=S.select_comparisons(space,{s['record']['inchikey'] for s in initial},'donor_rings')
    assert len(donors)==6
    for core in ('viologen','pmdi','ndi'):
        old=next(s for s in initial if s['record']['core']==core and s['contrast']=='ring_F')
        members=[s for s in donors if s['record']['core']==core]
        assert {s['contrast'] for s in members}=={'ring_OMe','ring_Me'}
        assert all(s['recipe']['template_atom']==old['recipe']['template_atom'] and
                   s['control']['inchikey']==old['control']['inchikey'] for s in members)
