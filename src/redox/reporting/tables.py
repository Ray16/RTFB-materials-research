"""Rebuild and independently audit tables in a new directory using existing records.

Each stage runs in a fresh interpreter. Legacy modules can use RESULTS globals
without redirecting another review or contaminating the caller's modules.
This workflow performs no quantum calculations and never promotes its output.

    python -m redox.reporting.tables --out results/reviews/my_review
"""
from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
import hashlib
import importlib
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

from redox.core.common import ROOT, RESULTS, CALCS
from redox.core.io import atomic_write_json
from redox.core.hashing import sha256_file
from redox.execution.environment import worker_environment


@dataclass(frozen=True)
class TableStage:
    name: str
    target: str
    method: str = "main"
    paths: tuple[tuple[str, str], ...] = (("RESULTS", "."),)
    arguments: tuple[tuple[str, str], ...] = ()
    dependencies: tuple[str, ...] = ()


STAGES = (
    TableStage("potentials", "redox.properties.potentials"),
    TableStage("orop_benchmark", "redox.quality.orop", "aggregate_all"),
    TableStage("reorganization", "redox.properties.reorg", "aggregate"),
    TableStage("integrity", "redox.properties.integrity"),
    TableStage("stability", "redox.properties.stability"),
    TableStage("stability_validation", "redox.quality.stability", dependencies=("stability",)),
    TableStage("capacity", "redox.properties.capacity_and_proxies"),
    TableStage("outer_sphere", "redox.reporting.outer_sphere", paths=(),
               arguments=(("output", "lambda_outer.csv"),)),
    TableStage("benchmark", "redox.reporting.benchmark", paths=(),
               arguments=(("results", "."), ("output", "validation/benchmark_mecn.csv")),
               dependencies=("potentials",)),
    TableStage("scorecard", "redox.screening.scorecard", "build",
               dependencies=("potentials", "orop_benchmark", "reorganization", "integrity",
                             "stability_validation", "capacity", "outer_sphere", "benchmark")),
    TableStage("pareto", "redox.screening.pareto", dependencies=("scorecard",)),
    TableStage("workbook", "redox.reporting.workbook", paths=(),
               arguments=(("results", "."), ("xlsx_output", "Candidates_computed.xlsx"),
                          ("csv_output", "starting_candidates_summary.csv")),
               dependencies=("pareto",)),
    TableStage("audit", "scripts/validation/audit/recompute_axes.py",
               paths=(("RES", "."), ("OUT", "validation")), dependencies=("workbook",)),
)
BENCHMARK_INPUTS = ("orop_benchmark.csv", "d3tales_reorg_validation_prod/comparison.csv")

# ---------------------------------------------------------------------------------------------
# Incremental (opt-in) support.
#
# Per-candidate, id-keyed stages that build their table purely by READING existing records
# (no quantum work happens here). Each one iterates the pool through a single helper imported
# from redox.core.common (`iterator`); in a --subset run that name is patched to a filtered
# view so the stage only emits rows for the requested candidates. Those subset rows are then
# merged with the previous review's cached rows. `order` names the canonical row ordering a
# FULL run produces, so a merged table is byte-identical to a from-scratch rebuild.
#
# Everything NOT listed here (orop_benchmark, benchmark, capacity, outer_sphere,
# stability_validation, scorecard, pareto, workbook, audit) is pool-level or cheap and is
# always rebuilt in full over the current tree -- incremental mode never risks their numbers.
INCREMENTAL_STAGES = {
    "potentials": dict(outputs=("redox_potentials.csv", "multielectron_stability.csv"),
                       iterator="read_manifest", id_column="id", order="manifest"),
    "reorganization": dict(outputs=("reorganization.csv",),
                           iterator="group_ids", id_column="id", order="group_ids"),
    "integrity": dict(outputs=("state_integrity.csv",),
                      iterator="group_ids", id_column="id", order="group_ids"),
    "stability": dict(outputs=("stability_disproportionation.csv",),
                      iterator="group_ids", id_column="id", order="group_ids"),
}

# Reference / benchmark namespaces that live under calcs/ but are NOT pool candidates. A change
# to any of them invalidates every referenced row, so they count as a global change.
RESERVED_IDS = frozenset({"ferrocene"})


def candidate_of(path: str) -> str | None:
    """Candidate id a fingerprinted input belongs to, or None when the input is global.

    Per-candidate inputs live under ``calcs/dft/<id>/...`` or ``calcs/uma/<id>/...`` (including
    a record's referenced gradients and the reorg cross-point caches). Everything else -- the
    shared scientific code, config, certificates, the manifest, the OROP/D3TaLES benchmark
    evidence and the ferrocene reference -- is global: a change there can move any row, so it
    forces a full rebuild rather than a per-candidate skip.
    """
    parts = Path(path).parts
    for i in range(len(parts) - 2):
        if parts[i] == "calcs" and parts[i + 1] in ("dft", "uma"):
            cid = parts[i + 2]
            if cid in RESERVED_IDS or cid.endswith(".bak") or ".bak" in cid:
                return None
            return cid
    return None


def classify_changes(previous: dict, current: dict) -> dict:
    """Diff two input inventories into per-candidate and global change sets.

    `previous`/`current` map input path -> sha256 (the input_inventory.json format). A path has
    changed if its hash differs or it was added/removed. Candidate ids with any changed input
    are recomputed; a changed global input means nothing can be reused.
    """
    keys = set(previous) | set(current)
    changed = [k for k in keys if previous.get(k) != current.get(k)]
    by_candidate: dict[str, list[str]] = {}
    global_changed: list[str] = []
    for key in changed:
        cid = candidate_of(key)
        if cid is None:
            global_changed.append(key)
        else:
            by_candidate.setdefault(cid, []).append(key)
    previous_ids = {cid for cid in (candidate_of(k) for k in previous) if cid is not None}
    changed_ids = sorted(cid for cid in by_candidate if cid in previous_ids)
    new_ids = sorted(cid for cid in by_candidate if cid not in previous_ids)
    return dict(changed_ids=changed_ids, new_ids=new_ids,
                recompute_ids=sorted(set(changed_ids) | set(new_ids)),
                global_changed=sorted(global_changed),
                global_changed_any=bool(global_changed),
                candidate_changed=bool(by_candidate))


def _read_csv(path: Path):
    with Path(path).open(newline="") as handle:
        reader = csv.DictReader(handle)
        return reader.fieldnames, list(reader)


def merge_csv(out_path: Path, base_path: Path | None, subset_paths, *, id_column: str,
              refreshed_ids, canonical_ids) -> int:
    """Write `out_path` = cached rows (for ids NOT refreshed) + freshly computed subset rows.

    Rows are read back as text and re-emitted, so the file is byte-identical to a full rebuild
    whenever the cached and fresh rows agree (they must, for an id whose inputs are unchanged).
    `refreshed_ids` are the candidates whose cached rows are DROPPED and replaced by subset rows;
    an id requested but absent from every subset file keeps its cached rows. `canonical_ids` is
    the full-run ordering of ids for this table.
    """
    refreshed = set(refreshed_ids)
    fieldnames = None
    rows_by_id: dict[str, list[dict]] = {}
    if base_path is not None and Path(base_path).is_file():
        fieldnames, base_rows = _read_csv(base_path)
        for row in base_rows:
            cid = row.get(id_column)
            if cid in refreshed:
                continue
            rows_by_id.setdefault(cid, []).append(row)
    for subset_path in subset_paths:
        if not Path(subset_path).is_file():
            continue
        subset_fields, subset_rows = _read_csv(subset_path)
        fieldnames = fieldnames or subset_fields
        for row in subset_rows:
            rows_by_id.setdefault(row.get(id_column), []).append(row)
    if fieldnames is None:
        raise FileNotFoundError(f"no cached or freshly computed rows to merge into {out_path}")
    canonical = list(canonical_ids)
    ordered = [cid for cid in canonical if cid in rows_by_id]
    leftover = sorted(cid for cid in rows_by_id if cid not in set(canonical))
    count = 0
    with Path(out_path).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for cid in [*ordered, *leftover]:
            for row in rows_by_id[cid]:
                writer.writerow(row)
                count += 1
    return count


def snapshot_inputs(source):
    """Fingerprint the global review's raw evidence and executing scientific code."""
    from redox.quality.provenance import referenced_files
    paths = {ROOT / "library/manifest.csv"}
    for directory in (CALCS, ROOT / "config", ROOT / "src/redox"):
        paths.update(path for path in directory.rglob("*") if path.is_file()
                     and path.suffix in (".json", ".xyz", ".npy", ".py")
                     and not path.name.startswith(("_opt_", "_seed_", "_tseed_"))
                     and "_superseded" not in path.parts and not any(".bak" in part for part in path.parts))
    # Per-record admission can depend on reviewed gradients/control artifacts
    # outside calcs/. The record's own bytes do not change if that evidence does.
    for record in tuple(paths):
        if record.suffix == ".json" and CALCS in record.parents:
            try:
                paths.update(referenced_files(json.loads(record.read_text()), ROOT))
            except (OSError, ValueError):
                # The raw record remains fingerprinted; readers refuse malformed
                # evidence. A later repair also changes that record's own hash.
                continue
    paths.update(ROOT / stage.target for stage in STAGES if stage.target.endswith(".py"))
    paths.update(Path(source) / name for name in BENCHMARK_INPUTS)
    paths.update((ROOT / "data/raw/validation/OROP").glob("*.csv"))
    paths.update((ROOT / "data/raw/validation/OROP/features").glob("*.csv"))
    paths.update((ROOT / "data/raw/validation/two_wave_mecn").glob("*.csv"))
    for directory in ("our_chunks", "our_chunks_v2"):
        paths.update((ROOT / "results/reorg_anchors" / directory).glob("chunk*_results.csv"))
    for certificate in (ROOT / "config").glob("native*.json"):
        paths.update(referenced_files(json.loads(certificate.read_text()), ROOT))
    model_review=ROOT/'config/model_validation_limits.json'
    if model_review.is_file():
        paths.update(referenced_files(json.loads(model_review.read_text()),ROOT))
    partition_certificate=ROOT/'config/projected_thermal_validation.json'
    if partition_certificate.is_file():
        paths.update(referenced_files(json.loads(partition_certificate.read_text()),ROOT))
    paths.update((ROOT/'scripts/validation/audit').glob('*.py'))
    for directory in (ROOT/'calcs').glob('**/projected_thermal'):
        paths.update(p for p in directory.rglob('*') if p.is_file() and p.suffix!='.lock')
    return {str(path): sha256_file(path) for path in paths if path.is_file()}


def _load_stage(stage: TableStage):
    if not stage.target.endswith(".py"):
        return importlib.import_module(stage.target)
    spec = importlib.util.spec_from_file_location(f"redox_table_{stage.name}", ROOT / stage.target)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _patch_subset(module, stage: TableStage, subset) -> None:
    """Restrict a per-candidate stage to `subset` ids by filtering the pool iterator it imported.

    Only the stage's OWN top-level iteration is narrowed; `read_result` and friends still consult
    the full manifest for each requested id (and the ferrocene reference), so every emitted row is
    byte-identical to the same row in a full run -- it is simply one of fewer rows.
    """
    name = INCREMENTAL_STAGES[stage.name]["iterator"]
    original = getattr(module, name)
    keep = set(subset)

    def filtered(*args, **kwargs):
        result = original(*args, **kwargs)
        if name == "read_manifest":
            return [row for row in result if row.get("id") in keep]
        return [gid for gid in result if gid in keep]

    setattr(module, name, filtered)


def _execute_stage(stage: TableStage, destination: Path, subset=None):
    """Used only in a disposable child interpreter."""
    module = _load_stage(stage)
    for attribute, relative in stage.paths:
        setattr(module, attribute, destination / relative)
    if subset is not None:
        if stage.name not in INCREMENTAL_STAGES:
            raise ValueError(f"stage {stage.name} does not support --subset")
        _patch_subset(module, stage, subset)
    result = getattr(module, stage.method)(**{key: destination / value for key, value in stage.arguments})
    if stage.name == "audit":
        atomic_write_json(destination / "stages/audit_result.json",
                          {"arithmetic_discrepancies": int(result)})
        return int(bool(result))
    return 0


def _canonical_ids(order: str) -> list[str]:
    """Full-run ordering of candidate ids for a per-record table (cheap; reads only disk)."""
    from redox.core.common import group_ids, read_manifest
    if order == "manifest":
        seen: dict[str, None] = {}
        for row in read_manifest():
            seen.setdefault(row["id"], None)
        return list(seen)
    return sorted(group_ids())


def _discover_previous(destination: Path) -> Path | None:
    """Newest sibling review directory usable as an incremental cache (most recent first)."""
    parent = Path(destination).resolve().parent
    if not parent.is_dir():
        return None
    candidates = []
    for child in parent.iterdir():
        if not child.is_dir() or child.resolve() == Path(destination).resolve():
            continue
        status, inventory = child / "review_status.json", child / "input_inventory.json"
        if not (status.is_file() and inventory.is_file()):
            continue
        try:
            state = json.loads(status.read_text()).get("status")
        except (OSError, ValueError):
            continue
        if state in ("complete", "reused", "inputs_changed", "candidate_failures"):
            candidates.append((status.stat().st_mtime, child))
    return max(candidates)[1] if candidates else None


def _reuse_previous(destination: Path, previous: Path, report: dict) -> None:
    """Copy a prior review's tables verbatim (byte-identical) when nothing changed."""
    skip = {"review_status.json", "input_inventory.json", "stages", "_incremental"}
    for item in Path(previous).iterdir():
        if item.name in skip:
            continue
        target = destination / item.name
        if item.is_dir():
            shutil.copytree(item, target, dirs_exist_ok=True)
        else:
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(item, target)
    try:
        prior = json.loads((Path(previous) / "review_status.json").read_text())
    except (OSError, ValueError):
        prior = {}
    report.update(status="reused", reused_from=str(previous),
                  arithmetic_discrepancies=prior.get("arithmetic_discrepancies"),
                  scientific_accuracy_validated=bool(prior.get("scientific_accuracy_validated")))


def review_tables(destination: Path, *, source: Path = RESULTS, workers: int = 1,
                  previous: Path | None = None, incremental: bool = False) -> dict:
    """Run CPU table stages, preserve existing reviews, and record every failure.

    With ``incremental=True`` the review is keyed on the per-candidate sha256 input inventory of
    a prior successful review (``previous``, or the newest usable sibling): unchanged candidates
    are reused byte-for-byte and only changed/new ones are recomputed. A single candidate's
    recompute failing no longer abandons the table -- it is recorded and the rest are written.
    The default call (``incremental=False``) is unchanged.
    """
    if type(workers) is not int or workers < 1:
        raise ValueError("workers must be a positive integer")
    cpu_sets = None
    if workers > 1:
        from redox.execution.resources import allocation_guard, cpu_allocation, host_capacity
        from redox.workflows.stages import validate_graph
        validate_graph(STAGES)
        if shutil.which("taskset") is None:
            raise RuntimeError("Parallel review requires taskset for disjoint CPU affinity")
        with allocation_guard():
            host_capacity(workers, minimum_ram_gib=2 * workers)
            cpu_sets = cpu_allocation(workers)
    destination = Path(destination).resolve()
    if incremental and previous is None:
        previous = _discover_previous(destination)
    previous = Path(previous).resolve() if previous is not None else None
    destination.mkdir(parents=True, exist_ok=False)
    logs = destination / "stages"
    logs.mkdir()
    report = dict(status="running", arithmetic_discrepancies=None,
                  scientific_accuracy_validated=False, quantum_jobs_submitted=0,
                  copied_benchmark_inputs=[], benchmark_input_sha256={}, stages=[],
                  workers=workers, cpu_affinity=cpu_sets, incremental=bool(incremental),
                  previous=str(previous) if previous else None, failed_candidates={})
    status_path = destination / "review_status.json"
    atomic_write_json(status_path, report)
    try:
        initial_inputs = snapshot_inputs(source)
        inventory = destination / "input_inventory.json"
        atomic_write_json(inventory, initial_inputs)
        report.update(input_inventory=str(inventory), inputs_current=None)
        for name in BENCHMARK_INPUTS:
            path = Path(source) / name
            if path.is_file():
                target = destination / name
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(path, target)
                report["copied_benchmark_inputs"].append(name)
                report["benchmark_input_sha256"][name] = sha256_file(target)

        # --- incremental decision ----------------------------------------------------------
        mode, recompute_ids = "full", []
        prior_inventory = (previous / "input_inventory.json") if previous else None
        if incremental and prior_inventory and prior_inventory.is_file():
            try:
                change = classify_changes(json.loads(prior_inventory.read_text()), initial_inputs)
            except (OSError, ValueError):
                change = None
            if change is not None:
                report["change_summary"] = change
                if not change["candidate_changed"] and not change["global_changed_any"]:
                    mode = "reuse"
                elif change["global_changed_any"]:
                    mode = "full"
                    report["incremental_fallback"] = "global_inputs_changed"
                else:
                    mode, recompute_ids = "partial", change["recompute_ids"]
        report["incremental_mode"] = mode

        if mode == "reuse":
            _reuse_previous(destination, previous, report)
            final_inputs = snapshot_inputs(source)
            changes = sorted(path for path in set(initial_inputs) | set(final_inputs)
                             if initial_inputs.get(path) != final_inputs.get(path))
            report.update(inputs_current=not changes, changed_inputs=changes)
            return report

        env = worker_environment(ROOT)
        entries = {}

        def started(stage):
            log = logs / f"{stage.name}.log"
            entry = dict(stage=stage.name, log=str(log), status="running", started=time.time())
            entries[stage.name] = entry
            report["stages"].append(entry)
            atomic_write_json(status_path, report)
            print(f"[tables] {stage.name}: {log}", flush=True)

        def _child(stage, slot, out_dir, subset_file=None):
            command = [sys.executable, "-m", "redox.reporting.tables", "--out", str(out_dir),
                       "--stage", stage.name]
            if subset_file is not None:
                command += ["--subset", str(subset_file)]
            if cpu_sets is not None:
                command = ["taskset", "--cpu-list", str(cpu_sets[slot]), *command]
            with (logs / f"{stage.name}.log").open("a") as stream:
                proc = subprocess.run(command, cwd=ROOT, env=env,
                                      stdout=stream, stderr=subprocess.STDOUT)
            return proc.returncode

        def _incremental_ready(stage):
            config = INCREMENTAL_STAGES[stage.name]
            return previous is not None and all((previous / name).is_file()
                                                for name in config["outputs"])

        def _merge_outputs(stage, subset_dirs, refreshed_ids):
            config = INCREMENTAL_STAGES[stage.name]
            canonical = _canonical_ids(config["order"])
            for name in config["outputs"]:
                merge_csv(destination / name, previous / name,
                          [Path(sub) / name for sub in subset_dirs],
                          id_column=config["id_column"],
                          refreshed_ids=refreshed_ids, canonical_ids=canonical)

        def run_incremental(stage, slot):
            base = destination / "_incremental" / stage.name
            base.mkdir(parents=True, exist_ok=True)
            subset_file = base / "subset.json"
            atomic_write_json(subset_file, recompute_ids)
            if _child(stage, slot, base, subset_file) == 0:
                _merge_outputs(stage, [base], recompute_ids)
                return 0
            # Isolate: a single bad candidate must not forfeit the others' rows.
            good, failed = [], []
            for cid in recompute_ids:
                one = base / f"id_{cid}"
                one.mkdir(parents=True, exist_ok=True)
                single = one / "subset.json"
                atomic_write_json(single, [cid])
                (good if _child(stage, slot, one, single) == 0 else failed).append((cid, one))
            good_ids = [cid for cid, _ in good]
            failed_ids = sorted(cid for cid, _ in failed)
            _merge_outputs(stage, [one for _, one in good], good_ids)
            if failed_ids:
                report["failed_candidates"].setdefault(stage.name, []).extend(failed_ids)
            return 0

        def run(stage, slot):
            if mode == "partial" and stage.name in INCREMENTAL_STAGES and _incremental_ready(stage):
                return run_incremental(stage, slot)
            return _child(stage, slot, destination)

        def finished(stage, returncode):
            entry = entries[stage.name]
            entry.update(returncode=returncode, status="failed" if returncode else "complete",
                         elapsed_seconds=time.time() - entry["started"])
            if stage.name == "audit" and (logs / "audit_result.json").exists():
                report.update(json.loads((logs / "audit_result.json").read_text()))
            if returncode:
                report.update(status="failed", failed_stage=stage.name)
            atomic_write_json(status_path, report)

        if workers == 1:
            # Keep exact historical stage ordering and support isolated stage tests.
            for stage in STAGES:
                started(stage)
                code = run(stage, 0)
                finished(stage, code)
                if code:
                    break
            else:
                report["status"] = "complete"
        else:
            from redox.workflows.stages import execute_graph
            execution = execute_graph(STAGES, run, workers=workers,
                                      started=started, finished=finished)
            report["stage_execution"] = execution
            if not execution["failed"]:
                report["status"] = "complete"
        if report["status"] == "complete" and report["failed_candidates"]:
            report["status"] = "candidate_failures"
        # A failed audit can coincide with a live calculation publishing new
        # input. Preserve both facts: do not hide its discrepancies, and do not
        # leave input freshness unknown merely because a stage returned nonzero.
        final_inputs = snapshot_inputs(source)
        changes = sorted(path for path in set(initial_inputs) | set(final_inputs)
                         if initial_inputs.get(path) != final_inputs.get(path))
        report.update(inputs_current=not changes, changed_inputs=changes)
        # Concurrent input drift is a stable-completion hazard for a from-scratch review, but in
        # incremental mode the next run simply recomputes the drifted candidate cheaply, so it is
        # recorded (not promoted to a failure) and never discards a review whose stages succeeded.
        if changes and report["status"] == "complete" and not incremental:
            report["status"] = "inputs_changed"
    except BaseException as exc:
        report.update(status="failed", error=repr(exc))
        raise
    finally:
        atomic_write_json(status_path, report)
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=1, help="concurrent independent CPU table stages; default 1")
    parser.add_argument("--incremental", action="store_true",
                        help="reuse unchanged candidates from a prior review; recompute only changed/new ones")
    parser.add_argument("--previous", type=Path, default=None,
                        help="prior review directory to reuse from (default: newest usable sibling of --out)")
    parser.add_argument("--stage", choices=[stage.name for stage in STAGES], help=argparse.SUPPRESS)
    parser.add_argument("--subset", type=Path, default=None, help=argparse.SUPPRESS)
    args = parser.parse_args(argv)
    if args.stage:
        stage = next(stage for stage in STAGES if stage.name == args.stage)
        subset = json.loads(args.subset.read_text()) if args.subset else None
        return _execute_stage(stage, args.out.resolve(), subset=subset)
    try:
        report = review_tables(args.out, workers=args.workers,
                               previous=args.previous, incremental=args.incremental)
    except FileExistsError:
        parser.error("choose a new output directory; existing reviews are preserved")
    print(f"[tables] {report['status']}: {args.out.resolve() / 'review_status.json'}")
    if report.get("incremental"):
        # Nonzero only when a stage failed or a requested candidate could not be built -- benign
        # concurrent input drift no longer discards a review whose work succeeded.
        return int(bool(report.get("failed_stage")) or bool(report.get("failed_candidates")))
    return int(report["status"] != "complete")


if __name__ == "__main__":
    raise SystemExit(main())
