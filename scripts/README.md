# Script entry points

Reusable implementations belong in `src/redox/`. Many older files here still
contain campaign-specific logic. This index identifies where to start without
renaming paths referenced by logs, manifests or validation evidence.

| Directory | Purpose | Entry points |
| --- | --- | --- |
| `pipeline/` | Candidate orchestration, collection and publication | `run_motifs.py` plans by default; `review_tables.py` wraps `redox.workflows.tables`; `collect_verified_candidates.py` collects the existing proposal cohort; `update_verified_slides.py` builds slides from verified snapshots |
| `mining/` | Proposal generation, registration and coverage | `identify_candidates.py`, `prepare_next_batch.py`, `prepare_candidate_campaign.py`; inspect dated cohort and worker assumptions before using |
| `validation/` | Physical controls and independent audits | `audit/recompute_axes.py`, `audit_candidate_provenance.py`, native thermal/solution reorganization controls; many commands run QM |
| `analysis/` | Sensitivity and conformer diagnostics | Inspect their inputs and methods before running |
| `plotting/` | Figures from derived or verified data | Subdirectories match candidate, validation, reorganization and pipeline figures |
| `fleet/` | Compatibility launchers | DFT shards and `gated_job.sh` delegate to package resource workflows; old campaign queues stop with a migration message |
| `legacy/fleet/` | Preserved historical launcher source | Reference text only; use the current workflow for new computations |
| `polaris/` | Polaris environment and PBS entry points | Use the Polaris skill and scheduler workflow |

`finalize_after_dft.sh` writes canonical tables, workbook and figures from saved
calculations after tests/readiness checks. Its old `SKIP_REORG_CROSS=0` compute
switch is rejected before any stage. For a software review or arithmetic audit,
prefer `python -m redox.workflows.tables --out results/reviews/NAME`.

`check_env.py` imports installed libraries and checks binaries. It is an environment
smoke test, not a GPU reservation or proof that a calculation is physically valid.
`free_gpus.py` now delegates to `redox.workflows.resources` and the authoritative
system `gpu_reserve list`. It no longer uses per-project host or occupancy policy.

For a new local DFT/UMA batch, `python -m redox.workflows.launch --module dft --out
logs/runs/NAME` prints a plan. Add `--execute` to run all suitable free local GPUs.
For a single reviewed validation command use `redox.workflows.gpu --gpu INDEX
--out logs/runs/NAME -- COMMAND`. Both paths record preflight/ownership and require
the validated 32-GB-class devices for the full molecular workflow. GPU library
probes happen inside the reservation. New scientific results still need their
domain-specific validation checks after execution.

For new candidate evaluations, start with
[`pipeline/screen_grafting_routes.py`](pipeline/screen_grafting_routes.py) and the
[route evidence](../results/grafting/README.md) before launching full redox work.
The active polymer interaction pilot is documented at
[results/polymer_models/README.md](../results/polymer_models/README.md):
`prepare_pi_stacking.py` builds conformer ensembles, `stage_polymer_conformers.py`
publishes checked warm starts, `run_pi_stacking.py` manages the finite reserved
queue, and `collect_pi_stacking.py` collects stationary electronic records.
`validation/review_polymer_interactions.py` supplies the independent structural
gate; `plotting/polymer_models/update_interaction_figures.py` maintains its plots
and three-dimensional browser. Existing campaign plans must be inspected before
any restart; preparation/submission actions are deliberately distinct from
read-only result collection.
