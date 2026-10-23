# Exact PMDI site-2 baseline: qualified workspace continuation

[plan.json](plan.json) preserves the original 115-atom chemical graph, warm start,
numerical model, optimizer settings and physical checks. It follows a separately
[checked large-point gradient recovery](../../../../validation/smd_gradient_workspace_32GB_2026-10-23/README.md).
[launch.json](launch.json), [preflight.json](preflight.json),
[startup_verification.json](startup_verification.json) and
[workspace_energy_status.json](workspace_energy_status.json) distinguish the
actual running optimizer from a completed electronic result.

One task starts on a reserved lambda2 V100 with two disjoint CPU cores and all
native thread limits set to two. The measured conservative density-fitting
host-RAM requirement is 218.6 GB; the node had sufficient available RAM at launch.
Other occupied GPUs are unavailable. The worker uses the original optimization,
durable geometry checkpoints, active single-point energies, independent force,
topology and structural checks. It does not compute or qualify native G.

Producer: [polymer_workspace_energy.py](../../../../../scripts/pipeline/polymer_workspace_energy.py).
Raw output and geometry evaluations:
`calcs/polymer_models/2026-10-23_crest_refinement/baseline__design_pmdi_ethyl__occ_2_workspace_32GB/`.
This is a single MMFF-start feasibility branch alongside the wider CREST search,
not a conformer-converged ensemble or representative experimental loading.

The previous candidate-only observer was stopped with exact process-identity
verification to prevent duplicate baseline submission during the handoff. Its
history and all earlier H100 attempts are preserved. The other two occupied
baseline components and subsequent thermal qualification still need continuation;
the independent wider-search observer and all fifteen CPU searches remain active
at the recorded live-process check.

The candidate observer subsequently resumed with the post-pilot site-3 branch.
This site-2 execution remains a separate, explicitly recorded continuation; its
physical outputs must pass before native-G handoff.
