"""Incremental, decoupled review_tables: fingerprint-skip, byte-identical merge, failure isolation.

These exercise the per-candidate incremental path with SYNTHETIC input fingerprints and a mocked
stage subprocess. No quantum work, no DFT, no GPU -- only the scheduling/reporting logic.
"""
import csv
import io
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from redox.reporting import tables


# --- pure classification -------------------------------------------------------------------

def test_candidate_of_separates_pool_candidates_from_global_inputs():
    assert tables.candidate_of("/r/calcs/dft/viologen/neu/result.json") == "viologen"
    assert tables.candidate_of("/r/calcs/uma/quinone/ox1/result.json") == "quinone"
    assert tables.candidate_of("/r/calcs/dft/dmpz/reorg/ox_at_neu__h.json") == "dmpz"
    # reference, benchmark, shared code and config are all global (not a single candidate)
    assert tables.candidate_of("/r/calcs/dft/ferrocene/neu/result.json") is None
    assert tables.candidate_of("/r/calcs/orop/16/ox/result.json") is None
    assert tables.candidate_of("/r/src/redox/properties/potentials.py") is None
    assert tables.candidate_of("/r/config/electrolyte.py") is None
    assert tables.candidate_of("/r/calcs/dft/viologen.bak/neu/result.json") is None


def test_classify_changes_skips_unchanged_recomputes_changed_and_flags_global():
    previous = {"/r/calcs/dft/A/neu/result.json": "1",
                "/r/calcs/dft/B/neu/result.json": "1",
                "/r/src/redox/x.py": "c"}
    current = {"/r/calcs/dft/A/neu/result.json": "1",          # unchanged -> skip
               "/r/calcs/dft/B/neu/result.json": "2",          # changed   -> recompute
               "/r/calcs/dft/C/neu/result.json": "9",          # new       -> recompute
               "/r/src/redox/x.py": "c"}
    change = tables.classify_changes(previous, current)
    assert change["changed_ids"] == ["B"]
    assert change["new_ids"] == ["C"]
    assert change["recompute_ids"] == ["B", "C"]
    assert change["global_changed_any"] is False


def test_classify_changes_detects_global_change():
    previous = {"/r/src/redox/properties/potentials.py": "v1",
                "/r/calcs/dft/A/neu/result.json": "1"}
    current = {"/r/src/redox/properties/potentials.py": "v2",   # scientific code changed
               "/r/calcs/dft/A/neu/result.json": "1"}
    change = tables.classify_changes(previous, current)
    assert change["global_changed_any"] is True
    assert change["recompute_ids"] == []


# --- subset injection ----------------------------------------------------------------------

def test_patch_subset_narrows_group_ids_and_manifest_iterators():
    reorg = SimpleNamespace(group_ids=lambda root=None: ["A", "B", "C"])
    tables._patch_subset(reorg, _reorg_stage(), {"B", "C"})
    assert reorg.group_ids() == ["B", "C"]              # pool loop sees only requested ids

    potentials_stage = next(s for s in tables.STAGES if s.name == "potentials")
    manifest = [dict(id="A", state="neu"), dict(id="B", state="neu"), dict(id="B", state="ox")]
    pot = SimpleNamespace(read_manifest=lambda: manifest)
    tables._patch_subset(pot, potentials_stage, {"B"})
    # every manifest row for a requested id is kept (both of B's states); others dropped
    assert pot.read_manifest() == [dict(id="B", state="neu"), dict(id="B", state="ox")]


# --- byte-identical merge ------------------------------------------------------------------

def _write_rows(path, fieldnames, rows):
    with Path(path).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def test_merge_is_byte_identical_to_a_full_rebuild(tmp_path):
    fields = ["id", "couple", "lambda_i_eV"]
    # Canonical full-pool order; A carries two couples (multi-row per id must be preserved).
    full_rows = [dict(id="A", couple="ox1->neu", lambda_i_eV=0.30),
                 dict(id="A", couple="neu->red1", lambda_i_eV=0.31),
                 dict(id="B", couple="ox->neu", lambda_i_eV=0.4261),
                 dict(id="C", couple="ox->neu", lambda_i_eV=0.50)]
    full = tmp_path / "full.csv"
    _write_rows(full, fields, full_rows)

    cached = tmp_path / "cached.csv"            # previous review had A and B
    _write_rows(cached, fields, [r for r in full_rows if r["id"] in ("A", "B")])
    subset = tmp_path / "subset.csv"            # recompute produced fresh B (changed) and C (new)
    _write_rows(subset, fields, [r for r in full_rows if r["id"] in ("B", "C")])

    out = tmp_path / "merged.csv"
    tables.merge_csv(out, cached, [subset], id_column="id",
                     refreshed_ids=["B", "C"], canonical_ids=["A", "B", "C"])
    assert out.read_bytes() == full.read_bytes()


# --- end-to-end incremental orchestration (mocked stage subprocess) -------------------------

def _reorg_stage():
    return next(stage for stage in tables.STAGES if stage.name == "reorganization")


def _seed_previous(tmp_path, inventory, rows):
    previous = tmp_path / "prev"
    previous.mkdir()
    (previous / "input_inventory.json").write_text(json.dumps(inventory))
    (previous / "review_status.json").write_text(json.dumps(
        {"status": "complete", "arithmetic_discrepancies": 0,
         "scientific_accuracy_validated": True}))
    _write_rows(previous / "reorganization.csv", ["id", "couple", "lambda_i_eV"], rows)
    (previous / "Candidates_computed.xlsx").write_text("binary-workbook")
    return previous


def test_reuse_copies_prior_tables_byte_for_byte_when_nothing_changed(tmp_path, monkeypatch):
    inventory = {"/r/calcs/dft/A/neu/result.json": "1", "/r/calcs/dft/B/neu/result.json": "1"}
    previous = _seed_previous(tmp_path, inventory,
                              [dict(id="A", couple="ox->neu", lambda_i_eV=0.30),
                               dict(id="B", couple="ox->neu", lambda_i_eV=0.40)])
    monkeypatch.setattr(tables, "snapshot_inputs", lambda source: dict(inventory))
    monkeypatch.setattr(tables.subprocess, "run",
                        lambda *a, **k: pytest.fail("reuse must run no stage"))

    report = tables.review_tables(tmp_path / "review", source=tmp_path / "src",
                                  previous=previous, incremental=True)
    assert report["status"] == "reused"
    assert report["arithmetic_discrepancies"] == 0
    out = tmp_path / "review"
    assert (out / "reorganization.csv").read_bytes() == (previous / "reorganization.csv").read_bytes()
    assert (out / "Candidates_computed.xlsx").read_text() == "binary-workbook"


def test_partial_recomputes_only_changed_candidates_and_merges(tmp_path, monkeypatch):
    previous_inv = {"/r/calcs/dft/A/neu/result.json": "1", "/r/calcs/dft/B/neu/result.json": "1"}
    previous = _seed_previous(tmp_path, previous_inv,
                              [dict(id="A", couple="ox->neu", lambda_i_eV=0.30),
                               dict(id="B", couple="ox->neu", lambda_i_eV=0.40)])
    current_inv = {"/r/calcs/dft/A/neu/result.json": "1",   # A unchanged -> cached row reused
                   "/r/calcs/dft/B/neu/result.json": "2",   # B changed   -> recompute
                   "/r/calcs/dft/C/neu/result.json": "9"}   # C new       -> recompute

    monkeypatch.setattr(tables, "snapshot_inputs", lambda source: dict(current_inv))
    monkeypatch.setattr(tables, "STAGES", (_reorg_stage(),))
    monkeypatch.setattr(tables, "_canonical_ids", lambda order: ["A", "B", "C"])
    subsets_seen = []

    def fake_run(command, **kwargs):
        out = Path(command[command.index("--out") + 1])
        subset = json.loads(Path(command[command.index("--subset") + 1]).read_text())
        subsets_seen.append(sorted(subset))
        _write_rows(out / "reorganization.csv", ["id", "couple", "lambda_i_eV"],
                    [dict(id=cid, couple="ox->neu", lambda_i_eV=0.99) for cid in subset])
        return SimpleNamespace(returncode=0)

    monkeypatch.setattr(tables.subprocess, "run", fake_run)
    report = tables.review_tables(tmp_path / "review", source=tmp_path / "src",
                                  previous=previous, incremental=True)

    assert report["incremental_mode"] == "partial"
    assert subsets_seen == [["B", "C"]]                 # only changed/new, in one batch
    merged = list(csv.DictReader((tmp_path / "review" / "reorganization.csv").open()))
    assert [r["id"] for r in merged] == ["A", "B", "C"]
    assert merged[0]["lambda_i_eV"] == "0.3"            # A reused byte-for-byte from cache
    assert merged[1]["lambda_i_eV"] == "0.99"           # B recomputed
    assert report["failed_candidates"] == {}


def test_one_candidates_failure_does_not_abandon_the_table(tmp_path, monkeypatch):
    previous_inv = {"/r/calcs/dft/A/neu/result.json": "1", "/r/calcs/dft/B/neu/result.json": "1"}
    previous = _seed_previous(tmp_path, previous_inv,
                              [dict(id="A", couple="ox->neu", lambda_i_eV=0.30),
                               dict(id="B", couple="ox->neu", lambda_i_eV=0.40)])
    current_inv = {"/r/calcs/dft/A/neu/result.json": "1",
                   "/r/calcs/dft/B/neu/result.json": "2",   # recompute, will succeed
                   "/r/calcs/dft/C/neu/result.json": "9"}   # recompute, will FAIL to build

    monkeypatch.setattr(tables, "snapshot_inputs", lambda source: dict(current_inv))
    monkeypatch.setattr(tables, "STAGES", (_reorg_stage(),))
    monkeypatch.setattr(tables, "_canonical_ids", lambda order: ["A", "B", "C"])

    def fake_run(command, **kwargs):
        out = Path(command[command.index("--out") + 1])
        subset = sorted(json.loads(Path(command[command.index("--subset") + 1]).read_text()))
        if subset == ["B", "C"] or subset == ["C"]:       # batch and C-alone both fail
            return SimpleNamespace(returncode=1)
        _write_rows(out / "reorganization.csv", ["id", "couple", "lambda_i_eV"],
                    [dict(id=cid, couple="ox->neu", lambda_i_eV=0.99) for cid in subset])
        return SimpleNamespace(returncode=0)

    monkeypatch.setattr(tables.subprocess, "run", fake_run)
    report = tables.review_tables(tmp_path / "review", source=tmp_path / "src",
                                  previous=previous, incremental=True)

    # The rest of the table is written: A reused, B recomputed; only C is recorded as failed.
    merged = list(csv.DictReader((tmp_path / "review" / "reorganization.csv").open()))
    assert [r["id"] for r in merged] == ["A", "B"]
    assert report["failed_candidates"] == {"reorganization": ["C"]}
    assert report["status"] == "candidate_failures"
    # main()'s exit-code contract: nonzero because a requested candidate could not be built,
    # but NOT merely because of benign input drift.
    assert bool(report.get("failed_stage")) is False
