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
    current_benchmark = json.loads((ROOT / "results/validation/publication_benchmark/current.json").read_text())
    reviews = [Path(p) for p in current_benchmark["review_snapshot_sha256"] if Path(p).name == "summary.json"]
    if len(reviews) != 1:
        raise ValueError("Expected the unique independent review behind the current benchmark")
    sources = {
        "benchmark.json": "results/validation/publication_benchmark/current.json",
        "benchmark_panel.json": "results/validation/publication_benchmark/panel_v1.json",
        "benchmark_review.json": str(reviews[0].relative_to(ROOT)),
        "base_assisted.json": "results/grafting/base_assisted/current.json",
        "grafting_status.json": "results/grafting/status.json",
        "polymer_validation.json": "results/reviews/2026-10-22/spacing_native_complete/validation.json",
        "segment_decomposition.json": "results/polymer_models/interaction_decomposition.json",
        "scalability.json": "results/validation/scalability_2026-10-21/summary.json",
        "shared_restart.json": "results/validation/shared_restart_2026-10-21/summary.json",
        "route_evidence.json": "results/grafting/route_evidence.json",
        "exact_routes.json": "results/validation/exact_routes_2026-10-22/summary.json",
        "acid_base.json": "results/validation/acid_base_external_2026-10-22/comparison.json",
        "cpu_concurrency.json": "results/validation/cpu_concurrency_2026-10-22/summary.json",
        "fresh_checkout.json": "results/validation/fresh_checkout_2026-10-22/summary.json",
        "reference_bias.json": "results/validation/reference_bias_2026-10-22/summary.json",
        "proton_transfer.json": "results/reviews/2026-10-22/proton_transfer/summary.json",
        "proton_transfer_sources.json": "results/validation/proton_transfer_2026-10-22/reference_snapshot.json",
        "quantum_restart.json": "results/validation/quantum_restart_2026-10-22/summary.json",
        "neutral_basis.json": "results/reviews/2026-10-22/neutral_basis_complete/summary.json",
        "exact_case_join.json": "results/reviews/2026-10-22/end_to_end_initial_join/summary.json",
        "cached_reproduction.json": "results/validation/cached_thermochemistry_2026-10-22/reconstruction.json",
        "cached_states_manifest.json": "examples/cached_thermochemistry/data/manifest.json",
        "saved_trimer_QC.json": "calcs/validation/saved_trimer_final_2026-10-22/status.json",
        "saved_trimer_fine_steps.json": "calcs/validation/saved_trimer_small_steps_2026-10-22/status.json",
        "benchmark_points.csv": "results/figures/validation/retrospective_redox_benchmark.csv",
        "finite_ensemble.json": "results/reviews/2026-10-23/baseline_native_complete/summary.json",
        "finite_ensemble_states.json": "results/reviews/2026-10-23/baseline_native_complete/current_snapshot.json",
        "pmdi_capacity.json": "results/validation/pmdi_occupancy_2026-10-23/capacity_scenarios.json",
        "pmdi_precursor_review.json": "results/reviews/2026-10-23/pmdi_precursor_access/summary.json",
        "all_axes_policy.json": "config/candidate_feasibility.json",
        "all_axes_summary.json": "results/selection/all_axes_2026-10-23_capacity_release/summary.json",
        "all_axes_decisions.json": "results/selection/all_axes_2026-10-23_capacity_release/decisions.json",
        "segment_relaxation.json": "results/reviews/2026-10-23/segment_deformation_audit_v2/summary.json",
        "basicity_reference_origin.json": "results/validation/dipea_reference_origin_2026-10-23_checked/summary.json",
    }
    plots = {
        "benchmark.pdf": "results/figures/validation/retrospective_redox_benchmark.pdf",
        "grafting.pdf": "results/grafting/base_assisted/figures/base_comparison.pdf",
        "construction.pdf": "results/polymer_models/figures/model_construction.pdf",
        "decomposition.pdf": "results/polymer_models/figures/fragment_decomposition.pdf",
        "conditional_grafting.pdf": "results/figures/polymer_models/conditional_grafting.pdf",
        "cpu_concurrency.pdf": "results/figures/validation/cpu_descriptor_concurrency.pdf",
        "reference_bias.pdf": "results/figures/validation/reference_bias.pdf",
        "spacing_construction.pdf": "results/figures/polymer_models/spacing_construction.pdf",
        "finite_ensemble.pdf": "results/figures/polymer_models/conditional_ensemble_baseline.pdf",
        "pmdi_capacity.pdf": "results/figures/polymer_models/pmdi_capacity_bounds.pdf",
        "pmdi_construction.pdf": "results/figures/polymer_models/pmdi_occupancy_construction_compact.pdf",
        "deformation_components.pdf": "results/figures/polymer_models/segment_deformation_components.pdf",
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
