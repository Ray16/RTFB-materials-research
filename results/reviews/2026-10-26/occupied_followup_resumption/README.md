# Occupied-chain follow-up resumption

[Startup snapshot](summary.json) · [Software tests](tests.xml) ·
[Native cohort](../../../polymer_models/campaigns/2026-10-26_occupied_site2_native_followup/batch.json) ·
[Checked occupied pilot](../occupied_native_pilot/independent_review.json)

The occupied site-2 pilot had already finished and passed independent native
review. Its earlier launch controller stopped before observing completion; the
independent review, not that stale dispatch status, is the source for the pilot.

Three previously declared electronic cohorts stopped when a combined multi-host
inventory query timed out. The revised CPU coordinators query permitted hosts
individually, retain unknown hosts as unknown, and reuse the existing bound-aware
preflight, atomic GPU reservation and physical workers. A separate native
coordinator manages 40 declared same-constitution states across two backbone
configurations. It reuses only the identical reviewed pilot; other states wait
for independently checked electronic records at the same geometries. No energy,
Hessian or thermochemistry is transferred between conformers.

At the recorded startup snapshot, all four CPU coordinators had live, verified
PID/start identities, and no new scientific GPU worker was running. Twelve old
submission intents remain ownership-unverified and are preserved; the controllers
do not resubmit those entries. Other entries may proceed after their geometry,
resource and reservation gates pass. A live coordinator or prepared entry is not a
finished calculation. The dispatcher and inventory records under each batch are
live; [summary.json](summary.json) is an immutable as-of snapshot.

The twelve held intents were subsequently checked and reopened; see the
[ownership reconciliation](../occupied_intent_reconciliation/README.md).

This is a phthalimide-control continuation in the alternative smooth-surface
MeCN model. It does not close the conformer sampling review comment, establish
the full polymer cycle, or change primary-model candidate values. Bulk worker
outputs and startup logs remain on shared storage at the paths in the snapshot.
