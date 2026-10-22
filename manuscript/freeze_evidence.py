#!/usr/bin/env python
"""Copy a deliberate manuscript evidence epoch without replacing earlier snapshots."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import time

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--name", required=True)
    args = parser.parse_args()
    if not args.name or Path(args.name).name != args.name:
        parser.error("Use a single snapshot directory name")
    evidence = ROOT / "manuscript/evidence" / args.name
    figures = ROOT / "results/figures/manuscript" / args.name
    if evidence.exists() or figures.exists():
        raise RuntimeError("An existing evidence epoch is immutable")
    sources = {
        "benchmark.json": "results/validation/publication_benchmark/current.json",
        "benchmark_panel.json": "results/validation/publication_benchmark/panel_v1.json",
        "benchmark_review.json": "results/reviews/2026-10-21/publication_benchmark_native_complete/summary.json",
        "base_assisted.json": "results/grafting/base_assisted/current.json",
        "grafting_status.json": "results/grafting/status.json",
        "polymer_validation.json": "results/reviews/2026-10-22/spacing_native_complete/validation.json",
        "segment_decomposition.json": "results/polymer_models/interaction_decomposition.json",
        "scalability.json": "results/validation/scalability_2026-10-21/summary.json",
        "shared_restart.json": "results/validation/shared_restart_2026-10-21/summary.json",
        "route_evidence.json": "results/grafting/route_evidence.json",
        "exact_routes.json": "results/validation/exact_routes_2026-10-22/summary.json",
        "acid_base.json": "results/validation/acid_base_external_2026-10-22/comparison.json",
        "benchmark_points.csv": "results/figures/validation/retrospective_redox_benchmark.csv",
    }
    plots = {
        "benchmark.pdf": "results/figures/validation/retrospective_redox_benchmark.pdf",
        "grafting.pdf": "results/grafting/base_assisted/figures/base_comparison.pdf",
        "construction.pdf": "results/polymer_models/figures/model_construction.pdf",
        "decomposition.pdf": "results/polymer_models/figures/fragment_decomposition.pdf",
        "conditional_grafting.pdf": "results/figures/polymer_models/conditional_grafting.pdf",
    }
    mapping = {evidence / name: ROOT / source for name, source in sources.items()}
    mapping.update({figures / name: ROOT / source for name, source in plots.items()})
    identities = {str(p.relative_to(ROOT)): sha(p) for p in mapping.values()}
    evidence.mkdir(); figures.mkdir()
    for destination, source in mapping.items():
        shutil.copyfile(source, destination)
    if any(sha(ROOT / p) != expected for p, expected in identities.items()):
        raise RuntimeError("A source changed while capturing the epoch; do not publish this snapshot")
    manifest = dict(kind="immutable_draft_evidence_snapshot", captured_UTC=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                    source_commit=subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
                    source_sha256=identities,
                    snapshot_sha256={str(p.relative_to(ROOT)): sha(p) for p in mapping},
                    scope="Working manuscript; source epochs retained; pending physical validation explicit")
    (evidence / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    pointer = dict(manifest=str((evidence / "manifest.json").relative_to(ROOT)), figures=str(figures.relative_to(ROOT)))
    (ROOT / "manuscript/evidence/current.json").write_text(json.dumps(pointer, indent=2) + "\n")
    print(json.dumps(pointer))


if __name__ == "__main__":
    main()
