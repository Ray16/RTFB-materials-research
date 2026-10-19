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

from redox.core.common import ROOT, RESULTS, CALCS
from redox.core.io import atomic_write_json


@dataclass(frozen=True)
class TableStage:
    name: str
    target: str
    method: str = "main"
    paths: tuple[tuple[str, str], ...] = (("RESULTS", "."),)


STAGES = (
    TableStage("potentials", "redox.properties.potentials"),
    TableStage("orop_benchmark", "redox.validation.orop", "aggregate_all"),
    TableStage("reorganization", "redox.properties.reorg", "aggregate"),
    TableStage("integrity", "redox.properties.integrity"),
    TableStage("stability", "redox.properties.stability"),
    TableStage("stability_validation", "redox.validation.stability"),
    TableStage("capacity", "redox.properties.capacity_and_proxies"),
    TableStage("outer_sphere", "scripts/pipeline/consolidate_lambda_outer.py",
               paths=(("OUT", "lambda_outer.csv"),)),
    TableStage("benchmark", "scripts/validation/benchmark/score_benchmark_mecn.py",
               paths=(("RESULTS", "."), ("OUT", "validation/benchmark_mecn.csv"))),
    TableStage("scorecard", "redox.screening.scorecard", "build"),
    TableStage("pareto", "redox.screening.pareto"),
    TableStage("workbook", "scripts/pipeline/fill_candidates_xlsx.py",
               paths=(("RESULTS", "."), ("XLSX_OUT", "Candidates_computed.xlsx"),
                      ("CSV_OUT", "starting_candidates_summary.csv"))),
    TableStage("audit", "scripts/validation/audit/recompute_axes.py",
               paths=(("RES", "."), ("OUT", "validation"))),
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
    paths.update(ROOT / stage.target for stage in STAGES if stage.target.endswith(".py"))
    paths.update(Path(source) / name for name in BENCHMARK_INPUTS)
    paths.update((ROOT / "data/raw/validation/OROP").glob("*.csv"))
    paths.update((ROOT / "data/raw/validation/OROP/features").glob("*.csv"))
    for directory in ("our_chunks", "our_chunks_v2"):
        paths.update((ROOT / "results/reorg_anchors" / directory).glob("chunk*_results.csv"))
    for certificate in (ROOT / "config").glob("native*.json"):
        paths.update(referenced_files(json.loads(certificate.read_text()), ROOT))
    partition_certificate=ROOT/'config/projected_thermal_validation.json'
    if partition_certificate.is_file():
        paths.update(referenced_files(json.loads(partition_certificate.read_text()),ROOT))
    paths.update((ROOT/'scripts/validation/audit').glob('*.py'))
    for directory in (ROOT/'calcs').glob('**/projected_thermal'):
        paths.update(p for p in directory.rglob('*') if p.is_file() and p.suffix!='.lock')
    return {str(path): hashlib.sha256(path.read_bytes()).hexdigest() for path in paths if path.is_file()}


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
    result = getattr(module, stage.method)()
    if stage.name == "audit":
        atomic_write_json(destination / "stages/audit_result.json",
                          {"arithmetic_discrepancies": int(result)})
        return int(bool(result))
    return 0


def review_tables(destination: Path, *, source: Path = RESULTS) -> dict:
    """Run CPU table stages, preserve existing reviews, and record every failure."""
    destination = Path(destination).resolve()
    destination.mkdir(parents=True, exist_ok=False)
    logs = destination / "stages"
    logs.mkdir()
    report = dict(status="running", arithmetic_discrepancies=None,
                  scientific_accuracy_validated=False, quantum_jobs_submitted=0,
                  copied_benchmark_inputs=[], benchmark_input_sha256={}, stages=[])
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
                report["benchmark_input_sha256"][name] = hashlib.sha256(target.read_bytes()).hexdigest()
        env = dict(os.environ)
        env["PYTHONPATH"] = os.pathsep.join(filter(None, (str(ROOT / "src"), env.get("PYTHONPATH"))))
        for key in ("OMP_NUM_THREADS", "OMP_THREAD_LIMIT", "MKL_NUM_THREADS",
                    "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
            env[key] = "1"
        for stage in STAGES:
            log = logs / f"{stage.name}.log"
            entry = dict(stage=stage.name, log=str(log), status="running")
            report["stages"].append(entry)
            atomic_write_json(status_path, report)
            print(f"[tables] {stage.name}: {log}", flush=True)
            with log.open("w") as stream:
                proc = subprocess.run(
                    [sys.executable, "-m", "redox.workflows.tables", "--out", str(destination),
                     "--stage", stage.name], cwd=ROOT, env=env,
                    stdout=stream, stderr=subprocess.STDOUT,
                )
            entry.update(returncode=proc.returncode, status="failed" if proc.returncode else "complete")
            if stage.name == "audit" and (logs / "audit_result.json").exists():
                report.update(json.loads((logs / "audit_result.json").read_text()))
            if proc.returncode:
                report.update(status="failed", failed_stage=stage.name)
                break
        else:
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
    parser.add_argument("--stage", choices=[stage.name for stage in STAGES], help=argparse.SUPPRESS)
    args = parser.parse_args(argv)
    if args.stage:
        stage = next(stage for stage in STAGES if stage.name == args.stage)
        return _execute_stage(stage, args.out.resolve())
    try:
        report = review_tables(args.out)
    except FileExistsError:
        parser.error("choose a new output directory; existing reviews are preserved")
    print(f"[tables] {report['status']}: {args.out.resolve() / 'review_status.json'}")
    return int(report["status"] != "complete")


if __name__ == "__main__":
    raise SystemExit(main())
