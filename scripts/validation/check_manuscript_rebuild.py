#!/usr/bin/env python
"""Rebuild committed paper artifacts in an isolated export without raw QM records.

Uses the current Python environment. Requires the figure dependencies and Graphviz;
also requires XeLaTeX/BibTeX unless --no-compile is supplied. The export, logs and JSON
report are retained in a new temporary directory (or a new --work-dir directory).
Byte equality checks figures, structures, SI, CSV tables and LaTeX/BibTeX sources.
PDF compilation is checked separately: PDF bytes can contain variable timestamps.
This verifies rendering from committed evidence, not regeneration of that evidence
from raw calculations, manuscript synchronization, or submission readiness.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import subprocess
import sys
import tarfile
import tempfile


ROOT = Path(__file__).resolve().parents[2]


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def artifacts(root):
    paths = [root / "manuscript/figures" / name for name in (
        "fig1_grafting.png", "fig2_pipeline.png", "fig3_qc_funnel.png",
        "fig4_validation_parity.png", "fig4b_residuals.png",
        "fig5a_anolyte_pareto.png", "fig5b_catholyte_pareto.png",
        "fig6_dmpz_design_rule.png", "fig6b_dmpz_catholyte_space.png",
        "fig7_cell_voltage.png", "toc_graphic.png")]
    paths += [root / name for name in (
        "manuscript/si_tables/S3_benchmark_endpoints.csv",
        "manuscript/si_tables/S7_scored_candidates.csv",
        "data/manuscript/pool_structures.csv", "manuscript/supporting_information.md",
        "manuscript/jcim/manuscript.tex", "manuscript/jcim/references.bib",
        "manuscript/jcim/supporting_information.tex")]
    return sorted(paths)


def compare_artifacts(root, expected):
    checks = {}
    for name, digest in expected.items():
        path = root / name
        actual = sha256(path) if path.is_file() else None
        checks[name] = dict(expected_sha256=digest, actual_sha256=actual, matches=actual == digest)
    return checks


def export_commit(commit, destination, archive):
    with archive.open("wb") as stream:
        subprocess.run(["git", "archive", "--format=tar", commit], cwd=ROOT,
                       stdout=stream, check=True)
    with tarfile.open(archive) as source:
        for member in source:
            path = destination / member.name
            if not path.resolve().is_relative_to(destination.resolve()):
                raise ValueError(f"Unsafe archive path: {member.name}")
            if member.isdir():
                path.mkdir(parents=True, exist_ok=True)
            elif member.isfile():
                path.parent.mkdir(parents=True, exist_ok=True)
                with source.extractfile(member) as stream:
                    path.write_bytes(stream.read())
            else:
                raise ValueError(f"Only regular tracked files are supported: {member.name}")
    archive.unlink()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ref", default="HEAD", help="Committed revision to verify (default HEAD)")
    parser.add_argument("--work-dir", type=Path, help="New directory for isolated export and logs")
    parser.add_argument("--no-compile", action="store_true", help="Check sources without building the PDF")
    args = parser.parse_args()
    commit = subprocess.check_output(
        ["git", "rev-parse", "--verify", "--end-of-options", args.ref + "^{commit}"],
        cwd=ROOT, text=True).strip()
    if args.work_dir is None:
        work = Path(tempfile.mkdtemp(prefix="redox-rebuild-"))
    else:
        work = args.work_dir.resolve()
        work.mkdir(parents=True, exist_ok=False)
    export = work / "source"
    export.mkdir()
    report = dict(commit=commit, python=sys.version, packages={}, commands=[],
                  scope="Rebuild from committed evidence; raw quantum records are not revalidated",
                  pdf_compiled=False, byte_identical=False, success=False)
    for name in ("numpy", "pandas", "matplotlib", "rdkit", "graphviz", "pillow"):
        try:
            report["packages"][name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            report["packages"][name] = None
    print(f"Rebuilding {commit[:12]} in {work}", flush=True)
    try:
        export_commit(commit, export, work / "source.tar")
        if (export / "calcs").exists() or (export / "library").exists():
            raise ValueError("Export unexpectedly contains raw calculation or library directories")
        expected = {str(p.relative_to(export)): sha256(p) for p in artifacts(export)}
        # Removing each output ensures a script cannot pass by leaving an old artifact intact.
        for name in expected:
            (export / name).unlink()
        (export / "manuscript/jcim/manuscript.pdf").unlink(missing_ok=True)
        (export / "manuscript/jcim/supporting_information.pdf").unlink(missing_ok=True)
        env = os.environ.copy()
        env.pop("REDOX_DATA_ROOT", None)
        env.pop("XTB_BIN", None)
        env.update(PYTHONPATH=str(export / "src"), OMP_NUM_THREADS="1", MKL_NUM_THREADS="1",
                   OPENBLAS_NUM_THREADS="1", MPLBACKEND="Agg", MPLCONFIGDIR=str(work / "mpl"))
        scripts = ["scripts/manuscript/pool_structures.py"]
        scripts += [str(p.relative_to(export)) for p in sorted(
            export.glob("scripts/plotting/manuscript/fig*.py"))]
        scripts += [
            "scripts/manuscript/build_si.py",
            "scripts/manuscript/build_si_pdf.py",
            "scripts/manuscript/build_latex.py",
        ]
        for i, script in enumerate(scripts):
            command = [sys.executable, script]
            if args.no_compile and script.endswith(("build_latex.py", "build_si_pdf.py")):
                command.append("--no-compile")
            log_path = work / f"{i:02d}_{Path(script).stem}.log"
            print(f"  {script}", flush=True)
            with log_path.open("w") as log:
                result = subprocess.run(command, cwd=export, env=env, stdout=log, stderr=subprocess.STDOUT)
            report["commands"].append(dict(command=command, returncode=result.returncode, log=log_path.name))
            if result.returncode:
                raise RuntimeError(f"{script} failed; see {log_path}")
        report["artifacts"] = compare_artifacts(export, expected)
        report["byte_identical"] = all(row["matches"] for row in report["artifacts"].values())
        if not report["byte_identical"]:
            mismatches = [name for name, row in report["artifacts"].items() if not row["matches"]]
            raise RuntimeError("Generated artifacts differ: " + ", ".join(mismatches))
        if not args.no_compile:
            for rel_pdf in ("manuscript/jcim/manuscript.pdf", "manuscript/jcim/supporting_information.pdf"):
                pdf = export / rel_pdf
                if not pdf.is_file() or not pdf.stat().st_size:
                    raise RuntimeError(f"Compilation produced no PDF: {rel_pdf}")
            report["pdf_compiled"] = True
        report["success"] = True
    except (OSError, ValueError, RuntimeError, subprocess.CalledProcessError) as exc:
        report["error"] = str(exc)
        print(str(exc), file=sys.stderr)
    finally:
        report_path = work / "report.json"
        report_path.write_text(json.dumps(report, indent=2) + "\n")
        print(f"Report: {report_path}", flush=True)
    if report["success"]:
        print(f"All {len(report['artifacts'])} artifacts are byte-identical. "
              f"PDF compiled: {report['pdf_compiled']}.")
    return 0 if report["success"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
