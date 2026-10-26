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
        "Fc_functional_review.json": "results/validation/fc_functional_sensitivity_2026-10-25/independent_review.json",
        "neighbor_cycle.json": "results/reviews/2026-10-25/pmdi_neighbor_cycle_checked/summary.json",
        "neighbor_native_review.json": "results/reviews/2026-10-25/pmdi_pair_native_residual_cuda_class_followup/component/summary.json",
        "neighbor_publication.json": "results/reviews/2026-10-25/pmdi_pair_native_published/publication.json",
        "neighbor_baseline.json": "results/reviews/2026-10-25/pmdi_pair_native_published/baseline_cycle_snapshot.json",
        "neighbor_figure.json": "results/figures/polymer_models/pmdi_neighbor_grafting.json",
        "neighbor_geometry_figure.json": "results/figures/polymer_models/pmdi_neighbor_geometry.json",
        "neighbor_analysis_source.py": "scripts/validation/analyze_neighbor_grafting.py",
        "neighbor_plot_source.py": "scripts/plotting/polymer_models/plot_neighbor_grafting.py",
        "neighbor_geometry_plot_source.py": "scripts/plotting/polymer_models/plot_neighbor_geometry.py",
        "neighbor_archive.json": "results/releases/neighbor_cycle/current.json",
        "neighbor_archive_relocation.json": "results/releases/neighbor_cycle/relocation.json",
        "neighbor_archive_manifest.json": "calcs/releases/2026-10-26_neighbor_cycle_native_replay/manifest.json",
        "neighbor_archive_tests.xml": "results/releases/neighbor_cycle/negative_tests.xml",
        "neighbor_replay_source.py": "examples/neighbor_cycle_replay/reproduce.py",
        "neighbor_archive_producer.py": "scripts/release/prepare_neighbor_addendum.py",
        "neighbor_relocation_source.py": "scripts/release/validate_neighbor_relocation.py",
        "control_addendum.json": "results/releases/controls/current.json",
        "control_archive_relocation.json": "results/releases/controls/relocation.json",
        "control_full_replay.json": "results/releases/controls/full_replay.json",
        "control_electronic_replay.json": "results/releases/controls/electronic_replay.json",
        "control_example_manifest.json": "examples/control_replay/data/manifest.json",
        "control_relocation_source.py": "scripts/release/validate_control_relocation.py",
        **{f"Fc_functional_{method}_{state}.json":
           f"calcs/validation/fc_functional_sensitivity_2026-10-25/{method}__{state}/electronic_diagnostic.json"
           for method in ('wb97mv', 'b3lyp', 'pbe0') for state in ('neu', 'ox')},
        "Fc_functional_source.py": "scripts/validation/fc_functional_sensitivity.py",
        "Fc_functional_review_source.py": "scripts/validation/review_fc_functional_sensitivity.py",
        "proton_transfer.json": "results/reviews/2026-10-22/proton_transfer/summary.json",
        "proton_transfer_sources.json": "results/validation/proton_transfer_2026-10-22/reference_snapshot.json",
        "quantum_restart.json": "results/validation/quantum_restart_2026-10-22/summary.json",
        "neutral_basis.json": "results/reviews/2026-10-22/neutral_basis_complete/summary.json",
        "exact_case_join.json": "results/reviews/2026-10-22/end_to_end_initial_join/summary.json",
        "pmdi_first_site_grafting.json": "results/reviews/2026-10-25/pmdi_first_site_grafting_checked/summary.json",
        "pmdi_first_site_producer.py": "scripts/validation/analyze_first_site_grafting.py",
        "pmdi_first_graft_figure.json": "results/figures/polymer_models/pmdi_first_graft.json",
        "pmdi_first_graft_plot_source.py": "scripts/plotting/polymer_models/plot_first_site_grafting.py",
        "CDS_contact.json": "results/validation/CDS_contact_2026-10-25/alternating_conf005.json",
        "CDS_smooth_component.json": "results/validation/smooth_CDS_2026-10-25/pilot/summary.json",
        "CDS_smooth_component_review.json": "results/validation/smooth_CDS_2026-10-25/pilot/independent_review.json",
        "CDS_smooth_full_gradients.json": "results/validation/smooth_CDS_2026-10-25/full_gradient_panel/summary.json",
        "CDS_smooth_quantum_CHO.json": "results/validation/smooth_quantum_2026-10-25/independent_review.json",
        "CDS_smooth_quantum_NCl.json": "results/validation/smooth_quantum_NCl_2026-10-25/independent_review.json",
        "CDS_smooth_qualification.json": "results/validation/smooth_model_qualification_2026-10-25/qualification.json",
        "CDS_smoothness_figure.json": "results/figures/validation/CDS_surface_smoothness.json",
        "CDS_smooth_model_source.py": "scripts/pipeline/smooth_cds_control.py",
        "CDS_smooth_control_source.py": "scripts/validation/qualify_smooth_quantum.py",
        "CDS_smooth_NCl_control_source.py": "scripts/validation/qualify_smooth_NCl.py",
        "CDS_smooth_review_source.py": "scripts/validation/review_smooth_quantum.py",
        "CDS_smoothness_plot_source.py": "scripts/plotting/validation/plot_CDS_surface_smoothness.py",
        "CDS_smooth_literature.json": "results/reviews/2026-10-25/smooth_surface_literature/reference.json",
        "pmdi_pair_electronic.json": "results/polymer_models/campaigns/2026-10-24_GPU_electronic/baseline__design_pmdi_ethyl__occ_2_3/energy/tasks/baseline__design_pmdi_ethyl__occ_2_3__pose0/status.json",
        "cached_reproduction.json": "results/validation/cached_thermochemistry_2026-10-22/reconstruction.json",
        "cached_states_manifest.json": "examples/cached_thermochemistry/data/manifest.json",
        "saved_trimer_QC.json": "calcs/validation/saved_trimer_final_2026-10-22/status.json",
        "saved_trimer_fine_steps.json": "calcs/validation/saved_trimer_small_steps_2026-10-22/status.json",
        "benchmark_points.csv": "results/figures/validation/retrospective_redox_benchmark.csv",
        "finite_ensemble.json": "results/reviews/2026-10-23/baseline_native_complete/summary.json",
        "finite_ensemble_states.json": "results/reviews/2026-10-23/baseline_native_complete/current_snapshot.json",
        "spacing_conformers.json": "results/figures/polymer_models/spacing_conformers.json",
        "pmdi_capacity.json": "results/validation/pmdi_occupancy_2026-10-23/capacity_scenarios.json",
        "pmdi_precursor_review.json": "results/reviews/2026-10-23/pmdi_precursor_access/summary.json",
        "all_axes_policy.json": "config/candidate_feasibility.json",
        "all_axes_summary.json": "results/selection/all_axes_2026-10-23_capacity_release/summary.json",
        "all_axes_decisions.json": "results/selection/all_axes_2026-10-23_capacity_release/decisions.json",
        "segment_relaxation.json": "results/reviews/2026-10-23/segment_deformation_audit_v2/summary.json",
        "basicity_reference_origin.json": "results/validation/dipea_reference_origin_2026-10-23_checked/summary.json",
        "review_package.json": "results/releases/primary_review/current.json",
        "review_package_manifest.json": "results/releases/primary_review/package_manifest.json",
        "review_package_data.json": "calcs/releases/2026-10-23_scientific_review_v5/data/manifest.json",
        "review_package_replay.json": "results/validation/release_package_2026-10-23/qualified_v5/summary.json",
        "review_package_arithmetic.json": "results/validation/release_package_2026-10-23/qualified_v5/arithmetic.json",
        "review_package_hessians.json": "results/validation/release_package_2026-10-23/qualified_v5/hessians.json",
        "TOC_artwork.json": "results/figures/publication/graphical_abstract_2026-10-23_v2.json",
        "TOC_artwork_check.json": "results/figures/publication/graphical_abstract_check.json",
        "thermal_cutoff_sensitivity.json": "results/validation/thermal_cutoff_sensitivity_2026-10-23/summary.json",
        "thermal_cutoff_review.json": "results/validation/thermal_cutoff_sensitivity_2026-10-23/independent_review.json",
        "thermal_cutoff_figure.json": "results/figures/validation/thermal_cutoff_sensitivity_v2.json",
        "review_package_draft.json": "results/releases/primary_review/draft_release.json",
        "quantum_environment.json": "results/validation/gpu_environment_restore_2026-10-23/checked/summary.json",
        "quantum_environment_states.json": "results/validation/gpu_environment_restore_2026-10-23/checked/states.json",
        "quantum_environment_compatibility.json": "results/validation/gpu_environment_restore_2026-10-23/compatibility_plan.json",
        "quantum_environment_review_source.py": "scripts/validation/publish_quantum_deployment.py",
        "quantum_environment_install_source.py": "scripts/release/restore_quantum_environment.py",
        "quantum_environment_certificate_source.py": "scripts/release/certificate_closure.py",
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
        "spacing_conformers.pdf": "results/figures/polymer_models/spacing_conformers.pdf",
        "pmdi_capacity.pdf": "results/figures/polymer_models/pmdi_capacity_bounds.pdf",
        "pmdi_first_graft.pdf": "results/figures/polymer_models/pmdi_first_graft.pdf",
        "CDS_surface_smoothness.pdf": "results/figures/validation/CDS_surface_smoothness.pdf",
        "pmdi_neighbor_grafting.pdf": "results/figures/polymer_models/pmdi_neighbor_grafting.pdf",
        "pmdi_neighbor_geometry.pdf": "results/figures/polymer_models/pmdi_neighbor_geometry.pdf",
        "pmdi_construction.pdf": "results/figures/polymer_models/pmdi_occupancy_construction_compact.pdf",
        "deformation_components.pdf": "results/figures/polymer_models/segment_deformation_components.pdf",
        "graphical_abstract.pdf": "results/figures/publication/graphical_abstract.pdf",
        "thermal_cutoff.pdf": "results/figures/validation/thermal_cutoff_sensitivity.pdf",
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
