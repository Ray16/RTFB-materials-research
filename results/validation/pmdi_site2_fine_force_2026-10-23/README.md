# Fine-step full-force diagnostic for the saved PMDI site-2 geometry

[plan.json](plan.json) predeclares a separate fixed-geometry diagnostic following
the preserved failed 0.0005/0.001-bohr force check. The
[CDS component audit](../pmdi_site2_CDS_surface_2026-10-23/README.md) locates the
sampled discrepancy in the legacy surface-energy term. This does not by itself
qualify the total force at smaller steps.

The new job evaluates full r2SCAN-D4/def2-SVP/SMD MeCN energies and gradients,
with the existing grid-response and adaptive conservative CDS derivative.
It checks 0.0001, 0.00005 and 0.00001 bohr, the original two sampled components
and the current two largest gradient components. Three complete baseline solves
test repeatability. Predeclared limits are 1e-6 Ha/bohr for analytic/FD error and
finest-pair agreement, 1e-11 Ha for baseline-energy spread and 1e-7 Ha/bohr for
repeated-gradient differences. No production threshold is relaxed.

The [launch](launch.json), [preflight](preflight.json) and
[actual startup check](startup_verification.json) verify a reserved V100 on
lambda2, two disjoint CPU cores, all six thread caps equal to two and one resident
process on the assigned GPU. The process has actually saved quantum evaluations;
[status.json](status.json) distinguishes its completed calculation from its actual startup.
Each completed SCF/gradient is checkpointed before proceeding.

The exact execution profile is supported by the checked 28 GiB spin controls.
Original failed outputs remain untouched. Raw quantum records are stored at
`calcs/validation/pmdi_site2_fine_force_2026-10-23/`; output/source hashes are
collected when the diagnostic finishes. No final electronic receipt, native
Hessian, G or candidate admission is made by this wrapper.

Use [polymer_gradient_profile.py](../../../scripts/validation/polymer_gradient_profile.py)
to prepare/run this kind of scoped diagnostic. After completion,
[review_polymer_gradient_profile.py](../../../scripts/validation/review_polymer_gradient_profile.py)
reconstructs its energy pairs, coordinates, repeatability and declared criteria
with standard Python. It refuses a still-running diagnostic. Passing this review
would establish only the stated derivative-convergence evidence, requiring
separate protocol review before production use.

The adjacent source searches are independent. Read-only process-start checks
verify one live CREST predictor on lambda2 and fourteen on lambda13; the
`live_searches_*.json` records retain process identities and artifact sizes.
They are not completed ensemble or DFT results.

The job completed all 15 SCF evaluations and passed [independent reconstruction](independent_review.json). Maximum analytic/FD error is 2.42e-7 Ha/bohr, finest-pair difference 2.14e-7 Ha/bohr, baseline energy spread 1.82e-12 Ha and repeated-gradient difference 2.90e-10 Ha/bohr. Peak device memory was 19,792 MiB. A [separate electronic qualification capsule](../../polymer_models/campaigns/2026-10-23_crest_refinement/baseline__design_pmdi_ethyl__occ_2_fine_force_qualified/README.md) explicitly reviews the combined evidence under a named four-step policy. It retains the original failure and the unchanged force tolerance. Native Hessian/G and the candidate cycle remain pending.

The native handoff exposed an orchestration error before any quantum work:
a preparation child inherited its observer's occupied single-CPU mask. The
stopped observer and source were preserved. Preparation children now recover the
initial permitted mask before the normal allocator selects their own free CPU.
The observer stays pinned to its own CPU. A real
[preparation integration check](preparation_affinity_integration.json) passes
with zero GPU jobs and an unregistered target; that directory is a metadata
control, not another scientific campaign. The observer was resumed with exactly
its original deadline and without stopping any scientific worker.

The actual [site-2 native pilot](../../polymer_models/campaigns/2026-10-23_crest_refinement/baseline__design_pmdi_ethyl__occ_2_native_xc_surface_budget/README.md)
subsequently started. [native_startup_verification.json](native_startup_verification.json)
confirms the live reserved process, disjoint CPU pair and inherited caps. Native
G remains pending. The two independent single-site components have checked E;
the adjacent two-graft E and all three native G values remain unresolved.
