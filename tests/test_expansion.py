"""Matched-comparison proposals and registration retain exact chemical identity."""
import copy
import hashlib
import json
import runpy

import pytest
from rdkit import Chem

from redox.build.contracts import validate_group
from redox.screening import expansion as E


def test_six_hypotheses_have_unique_explicit_grafts_and_complete_state_ladders():
    identities = set()
    for group in E.designs():
        mol = validate_group(E.SCAFFOLD, group)
        identities.add(Chem.MolToInchiKey(mol))
        assert group["provenance"]["exact_grafted_experiment"] is False
        assert group["provenance"]["production_qm"] is False
        assert len(group["states"]) == 3
        assert not any(key in group for key in ("E1", "E2", "lambda", "solubility", "experimental_potential"))
        if "cyclopropyl" in group["id"]:
            assert any(len(ring) == 3 for ring in mol.GetRingInfo().AtomRings())
        if "2fluoroethyl" in group["id"]:
            assert sum(atom.GetAtomicNum() == 9 for atom in mol.GetAtoms()) == 1
        if "ptz" not in group["id"]:
            assert len(mol.GetSubstructMatches(Chem.MolFromSmarts("[#6]=[#8]"))) == 4
    assert len(identities) == 6


@pytest.fixture
def registration(tmp_path, monkeypatch):
    run = tmp_path / "proposals"; run.mkdir()
    proposals = run / "proposals.json"
    proposals.write_text(json.dumps(E.designs()))
    (run / "preparation.json").write_text(json.dumps(dict(proposals_sha256=hashlib.sha256(proposals.read_bytes()).hexdigest())))
    config = tmp_path / "custom.py"
    config.write_text(f"# Existing researcher comments\nSCAFFOLD={E.SCAFFOLD!r}\nGROUPS=[]\n")
    monkeypatch.setattr(E, "registry_definitions", lambda: [])
    return run, config


def test_registration_is_append_only_preserves_original_and_is_idempotent(registration):
    run, config = registration
    original = config.read_bytes()
    result = E.register(run, config)
    assert len(result["registered"]) == 6
    assert config.read_bytes().startswith(original)
    archives = list((run / "config_before").glob("*.py"))
    assert len(archives) == 1 and archives[0].read_bytes() == original
    assert len(runpy.run_path(str(config))["GROUPS"]) == 6
    after = config.read_bytes()
    assert E.register(run, config)["registered"] == []
    assert config.read_bytes() == after


def test_an_existing_id_cannot_be_reassigned_to_a_new_structure(registration):
    run, config = registration
    altered = copy.deepcopy(E.designs()[0]); altered["frag"] = altered["frag"].replace("C5CC5", "C")
    config.write_text(f"SCAFFOLD={E.SCAFFOLD!r}\nGROUPS={[altered]!r}\n")
    before = config.read_bytes()
    with pytest.raises(ValueError, match="ID belongs to another exact structure"):
        E.register(run, config)
    assert config.read_bytes() == before


def test_changed_proposal_content_cannot_be_registered_under_an_old_preparation(registration):
    run, config = registration
    (run / "proposals.json").write_text("[]")
    before = config.read_bytes()
    with pytest.raises(RuntimeError, match="proposals changed"):
        E.register(run, config)
    assert config.read_bytes() == before
