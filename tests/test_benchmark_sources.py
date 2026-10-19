"""Experimental catholytes retain charge direction, identity and source limitations."""
import hashlib

import pytest
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

from redox.core.common import load_config
from redox.workflows import tables


@pytest.fixture
def benchmark():
    return load_config("benchmark_mecn")


def _phenothiazines(benchmark):
    return [b for b in benchmark.BENCHMARK
            if b["events"][0]["doi"] == "10.1021/acs.chemmater.8b04770"]


def test_oxidation_ladder_has_two_reduction_events_and_negative_signed_spacing(benchmark):
    b = next(b for b in _phenothiazines(benchmark)
             if b["name"] == "N-ethyl-3,7-dimethoxyphenothiazine")
    assert b["states"] == [("neu", 0, 1, 0), ("ox1", 1, 2, 1), ("ox2", 2, 1, 2)]
    assert [(e["event"], e["exp_V_vs_Fc"]) for e in b["events"] if e["kind"] == "E"] == [
        ("ox1->neu", .056), ("ox2->ox1", .655)]
    spacing = next(e for e in b["events"] if e["kind"] == "spacing")
    assert spacing["exp_V"] == pytest.approx(-.599)
    assert spacing["tier"] == "A"
    charges = {s[0]: s[1] for s in b["states"]}
    for e in b["events"]:
        if e["kind"] == "E":
            oxidized, reduced = e["event"].split("->")
            assert charges[oxidized] - charges[reduced] == 1


def test_source_conflicts_and_irreversibility_cannot_enter_strict_scoring(benchmark):
    parents = _phenothiazines(benchmark)
    assert len(parents) == 5
    strict = [(b, e) for b in parents for e in b["events"]
              if e["kind"] == "E" and e["tier"] == "A"]
    assert len(strict) == 5
    assert len({b["id"] for b, e in strict}) == 4
    conflicted = next(b for b in parents if b["name"].startswith("N-ethyl-3,7-bis"))
    assert {e["tier"] for e in conflicted["events"]} == {"B"}
    assert all("disagree" in e["tier_reason"] for e in conflicted["events"])
    excluded = [e for b in parents for e in b["events"] if e["tier"] == "excluded"]
    assert len(excluded) == 3
    assert all(e["event"] == "ox2->ox1" for e in excluded)


def test_second_oxidation_requires_an_explicit_charge_couple(benchmark):
    with pytest.raises(ValueError, match="Second oxidation lacks"):
        benchmark._layout(dict(name="unassigned oxidation", smiles_neutral="CC",
                               couple="+1/0", E1_vs_Fc_V=".1", E2_vs_Fc_V=".5"))


def test_one_wave_radical_oxidation_keeps_original_state_names_and_parity(benchmark):
    states, events = benchmark._layout(dict(name="TEMPO", smiles_neutral="[O]N1CCCCC1",
                                            couple="+1/0", E1_vs_Fc_V=".249", E2_vs_Fc_V=""))
    assert states == [("neu", 0, 2, 0), ("ox", 1, 1, 1)]
    assert events == [("+1/0", "ox->neu", .249, 1)]


def test_transcribed_phenothiazines_have_correct_formula_and_para_substitution(benchmark):
    formulas = {"N-ethylphenothiazine": "C14H13NS",
                "N-(2-(2-methoxyethoxy)ethyl)phenothiazine": "C17H19NO2S",
                "N-ethyl-3,7-dimethoxyphenothiazine": "C16H17NO2S",
                "3,7-dimethoxy-N-(2-(2-methoxyethoxy)ethyl)phenothiazine": "C19H23NO4S",
                "N-ethyl-3,7-bis(2-(2-methoxyethoxy)ethoxy)phenothiazine": "C24H33NO6S"}
    for b in _phenothiazines(benchmark):
        mol = Chem.MolFromSmiles(b["smiles"])
        assert rdMolDescriptors.CalcMolFormula(mol) == formulas[b["name"]]
        nitrogen = next(a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 7)
        aryl_ether_carbons = [n.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 8
                             for n in a.GetNeighbors() if n.GetIsAromatic()]
        assert len(aryl_ether_carbons) == (2 if "3,7" in b["name"] else 0)
        assert all(len(Chem.GetShortestPath(mol, nitrogen, c)) - 1 == 4
                   for c in aryl_ether_carbons)


def test_review_fingerprint_detects_changes_to_experimental_csv(tmp_path, monkeypatch):
    monkeypatch.setattr(tables, "ROOT", tmp_path)
    monkeypatch.setattr(tables, "CALCS", tmp_path / "calcs")
    source = tmp_path / "data/raw/validation/two_wave_mecn/phenothiazines_2019.csv"
    source.parent.mkdir(parents=True)
    source.write_text("name,E1_vs_Fc_V\nEPT,0.278\n")
    before = tables.snapshot_inputs(tmp_path / "results")
    assert before[str(source)] == hashlib.sha256(source.read_bytes()).hexdigest()
    source.write_text("name,E1_vs_Fc_V\nEPT,0.279\n")
    assert tables.snapshot_inputs(tmp_path / "results")[str(source)] != before[str(source)]
