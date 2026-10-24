# Exclusive-device two-graft checkpoint continuation

The preceding optimizer was stopped after an actual shared-GPU audit; the foreign
process was not disturbed. Its evaluated trajectory and 23 coupled energy/gradient
records remain intact. [checkpoint_review.json](checkpoint_review.json) identifies
the selected complete frame and checks exact graph/stereochemistry. The starting
force maximum is about 1.77e-4 Ha/bohr; this is an unqualified warm start, not an
optimized electronic result or free energy.

[plan.json](plan.json) retains the exact task, chemical model, validated execution
controls and all final gates. The optimizer state restarts in fresh output
storage. Prior same-graph H100 hybrid resource evidence is reused explicitly;
it is not presented as a new hybrid calculation or an exact-checkpoint resource
benchmark. The host allowance is doubled to cover potential screening changes.
No original calculation/status file or stale claim is removed.

The [finite residency observer](../../../../validation/gpu_exclusivity_2026-10-23/README.md)
is verified live before preparation, and watches actual future GPU residency.
The campaign is registered for an exclusively available H100, subject to shared
RAM/CPU checks. An occupied GPU cannot be used. Preparation is distinct from an
actual predictor start. Full optimization, active energies, structural/force
checks and native thermochemistry are still required. Candidate eligibility and
conformer convergence remain pending.

The subsequent actual start was interrupted by the residency guard before an accepted result: a delayed foreign CUDA attachment was confirmed. Only the owned process was stopped; no peer was signaled. The [second interruption record](../../../../validation/gpu_exclusivity_2026-10-23/second_residency_interruption.json) preserves the exact identities. This attempt is terminal and will not be automatically reused.
