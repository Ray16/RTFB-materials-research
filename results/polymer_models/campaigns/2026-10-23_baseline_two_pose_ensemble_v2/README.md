# Native free energies for both original spacing poses

This campaign preserves the original 12 checked electronic conformers and six
selected native free energies, then requests native thermochemistry for the other
six poses. It is a phthalimide method control, not a synthesis shortlist or a
substitute for the broader conformer/tacticity searches.

- [Frozen plan](plan.json): exact original tasks, input/producer hashes and resource bounds.
- [Current finite partitions](current.json): native counts, explicit missing G, duplicate sensitivity and conditional cycles.
- [Node preflight](lambda13_preflight.json), [actual start check](start_verification.json) and [worker record](thermal/workers/lambda13_gpu0_pilot.json).
- [Finite coordinator](coordinator.json): pilot-dependent continuation with a 24-hour window.
- Raw new states: `calcs/polymer_models/2026-10-23_baseline_two_pose_ensemble_v2/thermal/`.

The initial collector checked all six inherited native parents against their
active electronic records, exact saved geometries and the current native
thermochemistry protocol. All six two-pose RMSDs exceed 3 Å, so they are distinct
under the declared numerical duplicate criteria. The finite ensembles remain
blank until the missing G values pass their normal physical checks.

The first prepared directory, `2026-10-23_baseline_two_pose_ensemble`, is preserved
on shared storage as an unsubmitted preparation. A dispatch test found that a
resource pilot could otherwise choose an already completed energy-leading pose.
The corrected successor excludes completed thermal states from pilot selection;
its actual pilot is the missing `occ_2_4__pose0`. No QM job was launched from the
earlier preparation.

The completed panel now has all twelve native free energies and six finite
two-pose partitions. The [independent reconstruction](../../../reviews/2026-10-23/baseline_native_complete/summary.json)
also checks the balanced cycles and 25/50/100 cm⁻¹ rotor sensitivity. The old
coordinator stopped after a transient device-identity query timed out following
four successful worker starts. Scientific jobs continued and finished without
restarting. [Continuation observation](continuation_status.json) records that
separately from the original monitor failure. Broader sampling remains open.

See [the ensemble method and limits](../../../../docs/CONFORMER_ENSEMBLES.md).
