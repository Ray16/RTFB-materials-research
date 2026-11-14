"""Pinned calibration evidence: integrity, no-overwrite seeding, sigma provenance and the fail-loud guard."""
import json

import pytest

from redox.core.common import ROOT
from redox.screening import scorecard
from redox.workflows import reference


def test_bundle_matches_manifest_and_is_safe():
    manifest, members = reference.load()
    assert set(members) == set(manifest["files"])
    assert all(not name.startswith("/") and ".." not in name.split("/") for name in members)
    assert {row["id"] for row in manifest["manifest_rows"]} == set(manifest["ids"])
    assert "ferrocene" in manifest["ids"]


def test_seed_plan_refuses_to_overwrite(tmp_path):
    manifest, members = reference.load()
    name = sorted(members)[0]
    (tmp_path / name).parent.mkdir(parents=True)
    (tmp_path / name).write_bytes(members[name] + b"changed")
    present, missing, conflicts = reference.plan(manifest, members, root=tmp_path)
    assert conflicts == [name] and not present and len(missing) == len(members) - 1


def test_manifest_rows_must_match(tmp_path):
    manifest, _ = reference.load()
    row = dict(manifest["manifest_rows"][0])
    assert reference.merge_manifest_rows([row], root=tmp_path) == 1
    assert reference.merge_manifest_rows([row], root=tmp_path) == 0
    with pytest.raises(RuntimeError, match="differs"):
        reference.merge_manifest_rows([dict(row, smiles="C")], root=tmp_path)


def _evidence(redox="derived", disp="derived", lam="derived"):
    return {"sigma_sources": {"redox": redox, "disp": disp, "lambda": lam, "charge_class": {}}}


def test_prior_sigmas_are_refused_unless_allowed(monkeypatch):
    monkeypatch.delenv("REDOX_ALLOW_PRIOR_SIGMAS", raising=False)
    assert scorecard.require_derived_sigmas(_evidence()) == []
    with pytest.raises(RuntimeError, match="disp"):
        scorecard.require_derived_sigmas(_evidence(disp="prior"))
    monkeypatch.setenv("REDOX_ALLOW_PRIOR_SIGMAS", "1")
    assert scorecard.require_derived_sigmas(_evidence(redox="prior")) == ["redox"]


def test_missing_evidence_reports_prior_sources(monkeypatch, tmp_path):
    monkeypatch.setattr(scorecard, "RESULTS", tmp_path)
    sources = scorecard.derive_uncertainties([])["sigma_sources"]
    assert sources["redox"] == sources["disp"] == sources["lambda"] == "prior"


@pytest.mark.skipif(
    str(ROOT) != json.loads((reference.REFERENCE / f"{reference.DEFAULT}.manifest.json").read_text())["source_root"]
    or reference.plan(*reference.load())[1],
    reason="needs a seeded checkout at the records' origin (python -m redox seed-reference)")
def test_seeded_checkout_reproduces_validated_evidence():
    manifest, _ = reference.load()
    assert reference.verify(manifest) == []
