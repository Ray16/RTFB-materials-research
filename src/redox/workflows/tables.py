"""Rebuild and independently audit tables in a new directory using existing records.

Each stage runs in a fresh interpreter. Legacy modules can use RESULTS globals
without redirecting another review or contaminating the caller's modules.
This workflow performs no quantum calculations and never promotes its output.

    python -m redox.workflows.tables --out results/reviews/my_review
"""
from __future__ import annotations

import argparse
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
from redox.workflows.environment import worker_environment


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
    TableStage("orop_benchmark", "redox.validation.orop", "aggregate_all"),
    TableStage("reorganization", "redox.properties.reorg", "aggregate"),
    TableStage("integrity", "redox.properties.integrity"),
    TableStage("stability", "redox.properties.stability"),
    TableStage("stability_validation", "redox.validation.stability", dependencies=("stability",)),
    TableStage("capacity", "redox.properties.capacity_and_proxies"),
    TableStage("outer_sphere", "redox.reports.outer_sphere", paths=(),
               arguments=(("output", "lambda_outer.csv"),)),
    TableStage("benchmark", "redox.reports.benchmark", paths=(),
               arguments=(("results", "."), ("output", "validation/benchmark_mecn.csv")),
               dependencies=("potentials",)),
    TableStage("scorecard", "redox.screening.scorecard", "build",
               dependencies=("potentials", "orop_benchmark", "reorganization", "integrity",
                             "stability_validation", "capacity", "outer_sphere", "benchmark")),
    TableStage("pareto", "redox.screening.pareto", dependencies=("scorecard",)),
    TableStage("workbook", "redox.reports.workbook", paths=(),
               arguments=(("results", "."), ("xlsx_output", "Candidates_computed.xlsx"),
                          ("csv_output", "starting_candidates_summary.csv")),
               dependencies=("pareto",)),
    TableStage("audit", "scripts/validation/audit/recompute_axes.py",
               paths=(("RES", "."), ("OUT", "validation")), dependencies=("workbook",)),
)
BENCHMARK_INPUTS = ("orop_benchmark.csv", "d3tales_reorg_validation_prod/comparison.csv")


def snapshot_inputs(source):
    """Fingerprint the global review's raw evidence and executing scientific code."""
    from redox.validation.provenance import referenced_files
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


def _execute_stage(stage: TableStage, destination: Path):
    """Used only in a disposable child interpreter."""
    module = _load_stage(stage)
    for attribute, relative in stage.paths:
        setattr(module, attribute, destination / relative)
    result = getattr(module, stage.method)(**{key: destination / value for key, value in stage.arguments})
    if stage.name == "audit":
        atomic_write_json(destination / "stages/audit_result.json",
                          {"arithmetic_discrepancies": int(result)})
        return int(bool(result))
    return 0


def review_tables(destination: Path, *, source: Path = RESULTS, workers: int = 1) -> dict:
    """Run CPU table stages, preserve existing reviews, and record every failure."""
    if type(workers) is not int or workers < 1:
        raise ValueError("workers must be a positive integer")
    cpu_sets = None
    if workers > 1:
        from redox.workflows.resources import allocation_guard, cpu_allocation, host_capacity
        from redox.workflows.stages import validate_graph
        validate_graph(STAGES)
        if shutil.which("taskset") is None:
            raise RuntimeError("Parallel review requires taskset for disjoint CPU affinity")
        with allocation_guard():
            host_capacity(workers, minimum_ram_gib=2 * workers)
            cpu_sets = cpu_allocation(workers)
    destination = Path(destination).resolve()
    destination.mkdir(parents=True, exist_ok=False)
    logs = destination / "stages"
    logs.mkdir()
    report = dict(status="running", arithmetic_discrepancies=None,
                  scientific_accuracy_validated=False, quantum_jobs_submitted=0,
                  copied_benchmark_inputs=[], benchmark_input_sha256={}, stages=[],
                  workers=workers, cpu_affinity=cpu_sets)
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
        env = worker_environment(ROOT)
        entries = {}

        def started(stage):
            log = logs / f"{stage.name}.log"
            entry = dict(stage=stage.name, log=str(log), status="running", started=time.time())
            entries[stage.name] = entry
            report["stages"].append(entry)
            atomic_write_json(status_path, report)
            print(f"[tables] {stage.name}: {log}", flush=True)

        def run(stage, slot):
            command = [sys.executable, "-m", "redox.workflows.tables", "--out", str(destination),
                       "--stage", stage.name]
            if cpu_sets is not None:
                command = ["taskset", "--cpu-list", str(cpu_sets[slot]), *command]
            with (logs / f"{stage.name}.log").open("w") as stream:
                proc = subprocess.run(command, cwd=ROOT, env=env,
                                      stdout=stream, stderr=subprocess.STDOUT)
            return proc.returncode

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
        # A failed audit can coincide with a live calculation publishing new
        # input. Preserve both facts: do not hide its discrepancies, and do not
        # leave input freshness unknown merely because a stage returned nonzero.
        final_inputs = snapshot_inputs(source)
        changes = sorted(path for path in set(initial_inputs) | set(final_inputs)
                         if initial_inputs.get(path) != final_inputs.get(path))
        report.update(inputs_current=not changes, changed_inputs=changes)
        if changes and report["status"] == "complete":
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
    parser.add_argument("--stage", choices=[stage.name for stage in STAGES], help=argparse.SUPPRESS)
    args = parser.parse_args(argv)
    if args.stage:
        stage = next(stage for stage in STAGES if stage.name == args.stage)
        return _execute_stage(stage, args.out.resolve())
    try:
        report = review_tables(args.out, workers=args.workers)
    except FileExistsError:
        parser.error("choose a new output directory; existing reviews are preserved")
    print(f"[tables] {report['status']}: {args.out.resolve() / 'review_status.json'}")
    return int(report["status"] != "complete")


if __name__ == "__main__":
    raise SystemExit(main())
