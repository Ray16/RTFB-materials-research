# Migration

Source checkout: `/nfs/lambda_stor_01/homes/rzhu/0_redox` at
`e792621a0f140a627684ea411865d07d6e75e905`. This tree was built beside it.
The source files were not modified.

## Oracle comparison

All five portable replays match the pre-migration stdout byte for byte.

| Replay | Result |
|---|---|
| `examples/cached_thermochemistry/reproduce.py` | PMDI first reduction **-1.0757908637970104 V vs Fc**; DBU separated-ion reaction **-103.06882617264826 kJ/mol** |
| `examples/exact_viologen_replay/replay.py replay` | stdout sha256 `9ef6ba181c7568e438b7b04b94b79363770f5839ca448670cd5465ac44ae875e` |
| `examples/exact_viologen_corrected_replay/replay.py replay` | stdout sha256 `1e4f49e63807b85a75b0808e3a8f7ae7a530721bf8084f5be08814a5d6471c37` |
| `examples/finite_viologen_replay/replay.py` | stdout sha256 `f5d314956e55a531672059d4ddf79bf9977cf1440775be3c2dd59dae95ed928b` |
| `examples/occupied_sampling_replay/replay.py` | stdout sha256 `8a818b142c1684fd6a7ee5e20b235a44f5a9ecd738155f894a01b7b841a4a2f1` |

No scientific number was changed. Captured pre-migration stdout is in `baseline/`.

Curated tests: **744 passed, 1 skipped** (`pytest tests -q` in the `redox`
environment, 2026-11-16). The skip is the local xTB library guard in a test
that still ran its other cases. GPU parity tests were not added and no DFT
or Hessian job was launched.

## What the move did

Scientific functions were moved and their imports were rewritten. Equations,
constants, reduction order and admission rules were not edited, except for
path strings that point at this checkout.

- `build` became `structures`.
- `validation` became `quality`. Provenance readers `thermal_source` and
  `relocated_native` joined it.
- `reports` and `workflows/tables.py` became `reporting`.
- DFT, UMA, solvent, Hessian, response and checkpoint modules became `backends`.
  The four Hessian strategies and the two response roles remain separate modules.
- GFN2 diagnostics, vibrations, native thermochemistry, native repair and
  projected partitions became `thermochemistry`. `thermal.py` is
  `thermochemistry/xtb_diagnostic.py` so its name is not confused with native DFT.
- Resource discovery, GPU reservation, launch, chunk ownership and checkpoint
  transport became `execution`.
- Campaign and recovery modules stayed in `workflows`. Several of them locate
  launcher files as `src/redox/workflows/<name>.py`, and archived namespaces
  still refer to the historical archive paths. Those archive paths were not
  rewritten.

## Byte-identical files

| File | Why it cannot move or be edited |
|---|---|
| `src/redox/qm/qrrho.py` | Projected-thermal certificate admits only the xTB library-path transformation. |
| `src/redox/qm/xtb_kernel.py` | Sibling worker named by `qrrho.py` and by the certificate. |
| `src/redox/qm/xtb_runtime.py` | Import path is inside the certified `qrrho.py` bytes. Its own imports of native and projected thermochemistry were updated. |
| `src/redox/quality/sampling_snapshot.py` | Frozen occupied-sampling receipt stores this file's sha256 as `producer_sha256`. |
| `scripts/plotting/plot_occupied_sampling.py` | Frozen figure record stores this file's sha256. |
| Hash-locked example payloads | Replays refuse any content change. `examples/occupied_sampling_replay/replay.py` is not in that manifest and now imports `execution` and `quality`. |

Compatibility imports, used only so those frozen files still load:

- `redox.workflows.resources` → `redox.execution.resources`
- `redox.validation.sampling_receipts` → `redox.quality.sampling_receipts`
- `redox.validation.sampling_native_receipts` → `redox.quality.sampling_native_receipts`

## Left behind

Not copied: `manuscript/`, `slides/`, `validation/` campaign directories,
`logs/`, `library/`, `calcs/`, bulk `results/`, and
`scripts/report/.browsers/` plus `scripts/report/.check/`.

Pinned files required by a certificate or a replay were copied to the same
relative paths:

- `config/projected_thermal_validation.json` evidence under
  `results/validation/rigid_mode_thermal_2026-10-18/` and
  `results/reviews/2026-10-18/software_checks_isolated_kernel_full.log`
- `results/figures/slides/occupied_sampling_convergence.pdf`

`data/redox_route_sources/` was not copied. It is gitignored local source
material in the original checkout.

Tests that only opened those excluded trees were not kept:
`test_context_consequences.py`, `test_phenothiazine_route.py`,
`test_qualified_cases.py`, `test_sampling_batches.py`,
`test_screening_sensitivity.py`, `test_viologen_exact_route.py`, and the
manuscript, release-archive and dated-campaign tests that were never copied.
The occupied-sampling numerical replay still covers the batch comparison those
campaign fixtures exercised.

## Public commands

`redox` and `redox-uma` keep the same command names. Dispatch targets are the
new modules (`redox.backends.uma_current`, `redox.thermochemistry.native`,
`redox.execution.resources`, and the rest listed in `src/redox/cli.py`).

## Module map

| Old path | New path |
|---|---|
| `src/redox/build/__init__.py` | `src/redox/structures/__init__.py` |
| `src/redox/build/attachment_identity.py` | `src/redox/structures/attachment_identity.py` |
| `src/redox/build/candidates.py` | `src/redox/structures/candidates.py` |
| `src/redox/build/contracts.py` | `src/redox/structures/contracts.py` |
| `src/redox/build/geometry.py` | `src/redox/structures/geometry.py` |
| `src/redox/build/grafting.py` | `src/redox/structures/grafting.py` |
| `src/redox/build/ionic_polymer.py` | `src/redox/structures/ionic_polymer.py` |
| `src/redox/build/manifest.py` | `src/redox/structures/manifest.py` |
| `src/redox/build/oligomer_topology.py` | `src/redox/structures/oligomer_topology.py` |
| `src/redox/build/oligomers.py` | `src/redox/structures/oligomers.py` |
| `src/redox/build/polymer_cases.py` | `src/redox/structures/polymer_cases.py` |
| `src/redox/build/polymer_contacts.py` | `src/redox/structures/polymer_contacts.py` |
| `src/redox/build/standalone.py` | `src/redox/structures/standalone.py` |
| `src/redox/build/structures.py` | `src/redox/structures/structures.py` |
| `src/redox/build/validation.py` | `src/redox/structures/validation.py` |
| `src/redox/qm/cds_finite.py` | `src/redox/backends/cds_finite.py` |
| `src/redox/qm/compact_smd.py` | `src/redox/backends/compact_smd.py` |
| `src/redox/qm/dft.py` | `src/redox/backends/dft.py` |
| `src/redox/qm/gas_reference.py` | `src/redox/backends/gas_reference.py` |
| `src/redox/qm/hessian_memory.py` | `src/redox/backends/hessian_memory.py` |
| `src/redox/qm/hessian_restart.py` | `src/redox/backends/hessian_restart.py` |
| `src/redox/qm/hessian_resume.py` | `src/redox/backends/hessian_resume.py` |
| `src/redox/qm/native_thermal.py` | `src/redox/thermochemistry/native.py` |
| `src/redox/qm/native_thermal_repair.py` | `src/redox/thermochemistry/native_repair.py` |
| `src/redox/qm/nonequilibrium_pcm.py` | `src/redox/backends/nonequilibrium_pcm.py` |
| `src/redox/qm/nonequilibrium_smd.py` | `src/redox/backends/nonequilibrium_smd.py` |
| `src/redox/qm/occupied_checkpoints.py` | `src/redox/backends/occupied_checkpoints.py` |
| `src/redox/qm/pcm_bounded.py` | `src/redox/backends/pcm_bounded.py` |
| `src/redox/qm/projected_thermal.py` | `src/redox/thermochemistry/projected.py` |
| `src/redox/qm/relocated_native.py` | `src/redox/quality/relocated_native.py` |
| `src/redox/qm/response_diagnostic.py` | `src/redox/backends/response_diagnostic.py` |
| `src/redox/qm/response_resume.py` | `src/redox/backends/response_resume.py` |
| `src/redox/qm/sampling_native_checkpoints.py` | `src/redox/backends/sampling_native_checkpoints.py` |
| `src/redox/qm/saved_hessian.py` | `src/redox/backends/saved_hessian.py` |
| `src/redox/qm/scf_checkpoints.py` | `src/redox/backends/scf_checkpoints.py` |
| `src/redox/qm/screening_barrier.py` | `src/redox/backends/screening_barrier.py` |
| `src/redox/qm/sp.py` | `src/redox/backends/sp.py` |
| `src/redox/qm/stable_geometry_scf.py` | `src/redox/backends/stable_geometry_scf.py` |
| `src/redox/qm/thermal.py` | `src/redox/thermochemistry/xtb_diagnostic.py` |
| `src/redox/qm/thermal_source.py` | `src/redox/quality/thermal_source.py` |
| `src/redox/qm/torsion_scan.py` | `src/redox/backends/torsion_scan.py` |
| `src/redox/qm/uma.py` | `src/redox/backends/uma.py` |
| `src/redox/qm/uma_current.py` | `src/redox/backends/uma_current.py` |
| `src/redox/qm/uma_model_identity.py` | `src/redox/backends/uma_model_identity.py` |
| `src/redox/qm/uma_nitro.py` | `src/redox/backends/uma_nitro.py` |
| `src/redox/qm/uma_policy.py` | `src/redox/backends/uma_policy.py` |
| `src/redox/qm/vibrations.py` | `src/redox/thermochemistry/vibrations.py` |
| `src/redox/reports/__init__.py` | `src/redox/reporting/__init__.py` |
| `src/redox/reports/benchmark.py` | `src/redox/reporting/benchmark.py` |
| `src/redox/reports/outer_sphere.py` | `src/redox/reporting/outer_sphere.py` |
| `src/redox/reports/workbook.py` | `src/redox/reporting/workbook.py` |
| `src/redox/validation/__init__.py` | `src/redox/quality/__init__.py` |
| `src/redox/validation/archived_sources.py` | `src/redox/quality/archived_sources.py` |
| `src/redox/validation/benchmark_panel.py` | `src/redox/quality/benchmark_panel.py` |
| `src/redox/validation/cds_hessian_component.py` | `src/redox/quality/cds_hessian_component.py` |
| `src/redox/validation/charged_hessian_source.py` | `src/redox/quality/charged_hessian_source.py` |
| `src/redox/validation/component_sensitivity.py` | `src/redox/quality/component_sensitivity.py` |
| `src/redox/validation/contact_cohort.py` | `src/redox/quality/contact_cohort.py` |
| `src/redox/validation/context_consequences.py` | `src/redox/quality/context_consequences.py` |
| `src/redox/validation/evidence_epoch.py` | `src/redox/quality/evidence_epoch.py` |
| `src/redox/validation/execution_review.py` | `src/redox/quality/execution_review.py` |
| `src/redox/validation/failed_execution.py` | `src/redox/quality/failed_execution.py` |
| `src/redox/validation/ionic_preparation.py` | `src/redox/quality/ionic_preparation.py` |
| `src/redox/validation/model_limits.py` | `src/redox/quality/model_limits.py` |
| `src/redox/validation/molecular_native.py` | `src/redox/quality/molecular_native.py` |
| `src/redox/validation/namespace_sources.py` | `src/redox/quality/namespace_sources.py` |
| `src/redox/validation/native_component_sensitivity.py` | `src/redox/quality/native_component_sensitivity.py` |
| `src/redox/validation/native_conformer_batches.py` | `src/redox/quality/native_conformer_batches.py` |
| `src/redox/validation/orop.py` | `src/redox/quality/orop.py` |
| `src/redox/validation/portable_energy_parent.py` | `src/redox/quality/portable_energy_parent.py` |
| `src/redox/validation/portable_inputs.py` | `src/redox/quality/portable_inputs.py` |
| `src/redox/validation/primary_checkpoints.py` | `src/redox/quality/primary_checkpoints.py` |
| `src/redox/validation/protonation.py` | `src/redox/quality/protonation.py` |
| `src/redox/validation/provenance.py` | `src/redox/quality/provenance.py` |
| `src/redox/validation/readiness.py` | `src/redox/quality/readiness.py` |
| `src/redox/validation/reference_admission.py` | `src/redox/quality/reference_admission.py` |
| `src/redox/validation/reference_diagnostics.py` | `src/redox/quality/reference_diagnostics.py` |
| `src/redox/validation/release_gate.py` | `src/redox/quality/release_gate.py` |
| `src/redox/validation/sampling_batches.py` | `src/redox/quality/sampling_batches.py` |
| `src/redox/validation/sampling_native_receipts.py` | `src/redox/quality/sampling_native_receipts.py` |
| `src/redox/validation/sampling_publication.py` | `src/redox/quality/sampling_publication.py` |
| `src/redox/validation/sampling_receipts.py` | `src/redox/quality/sampling_receipts.py` |
| `src/redox/validation/sampling_snapshot.py` | `src/redox/quality/sampling_snapshot.py` |
| `src/redox/validation/scf_frame.py` | `src/redox/quality/scf_frame.py` |
| `src/redox/validation/screening_sensitivity.py` | `src/redox/quality/screening_sensitivity.py` |
| `src/redox/validation/solution_pair_review.py` | `src/redox/quality/solution_pair_review.py` |
| `src/redox/validation/solution_relaxation.py` | `src/redox/quality/solution_relaxation.py` |
| `src/redox/validation/source_inventory.py` | `src/redox/quality/source_inventory.py` |
| `src/redox/validation/spin_ensemble.py` | `src/redox/quality/spin_ensemble.py` |
| `src/redox/validation/spin_evidence.py` | `src/redox/quality/spin_evidence.py` |
| `src/redox/validation/spins.py` | `src/redox/quality/spins.py` |
| `src/redox/validation/stability.py` | `src/redox/quality/stability.py` |
| `src/redox/validation/stationarity_failure.py` | `src/redox/quality/stationarity_failure.py` |
| `src/redox/validation/uma_environment.py` | `src/redox/quality/uma_environment.py` |
| `src/redox/validation/uma_parity.py` | `src/redox/quality/uma_parity.py` |
| `src/redox/validation/workspace_index.py` | `src/redox/quality/workspace_index.py` |
| `src/redox/workflows/bounded_molecular.py` | `src/redox/execution/bounded_molecular.py` |
| `src/redox/workflows/checkpoint_io.py` | `src/redox/execution/checkpoint_io.py` |
| `src/redox/workflows/checkpoints.py` | `src/redox/execution/checkpoints.py` |
| `src/redox/workflows/chunked.py` | `src/redox/execution/chunked.py` |
| `src/redox/workflows/cpu_chunks.py` | `src/redox/execution/cpu_chunks.py` |
| `src/redox/workflows/environment.py` | `src/redox/execution/environment.py` |
| `src/redox/workflows/geometry_checkpoints.py` | `src/redox/execution/geometry_checkpoints.py` |
| `src/redox/workflows/gpu.py` | `src/redox/execution/gpu.py` |
| `src/redox/workflows/host_local_review.py` | `src/redox/execution/host_local_review.py` |
| `src/redox/workflows/launch.py` | `src/redox/execution/launch.py` |
| `src/redox/workflows/launch_current.py` | `src/redox/execution/launch_current.py` |
| `src/redox/workflows/owned_debugger.py` | `src/redox/execution/owned_debugger.py` |
| `src/redox/workflows/prelaunch_accounting.py` | `src/redox/execution/prelaunch_accounting.py` |
| `src/redox/workflows/prepared_namespace_queue.py` | `src/redox/execution/prepared_namespace_queue.py` |
| `src/redox/workflows/processes.py` | `src/redox/execution/processes.py` |
| `src/redox/workflows/publication.py` | `src/redox/execution/publication.py` |
| `src/redox/workflows/residency_refusal.py` | `src/redox/execution/residency_refusal.py` |
| `src/redox/workflows/resources.py` | `src/redox/execution/resources.py` |
| `src/redox/workflows/submission_refusal.py` | `src/redox/execution/submission_refusal.py` |
| `src/redox/workflows/tables.py` | `src/redox/reporting/tables.py` |
