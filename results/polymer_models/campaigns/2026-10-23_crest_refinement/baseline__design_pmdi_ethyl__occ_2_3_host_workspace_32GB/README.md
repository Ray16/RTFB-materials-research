# Full 139-atom active-energy reference point on V100

The [plan](plan.json) retains the exact starting geometry and quantum settings
of the completed large-device reference. The
[actual scientific start](startup_verification.json) is verified on lambda2
GPU 0, a 32 GiB V100, with disjoint CPUs 5 and 21 and native thread caps of two.
The [live footprint record](footprint_qualification.json) is running; no final
energy, optimized minimum or free energy is admitted from this start.

The full- and range-separated screened DF arrays require 163,190,172,552 and
68,435,226,384 bytes. The original inventory rule adds 20% plus 24 GiB, requiring
303,720,282,499 bytes of available host memory. The worker has a separate
283 GiB RSS bound and 32 GiB free-RAM floor. The established two-spin numerical
controls, 28 GiB workspace budget and 31 GiB measured GPU ceiling remain required.

The concurrently running geometry stage has a separately enforced 64 GiB host
bound, above its measured 1.5 GiB peak. Its existing optimizer, checkpoint and
scientific settings are untouched. The new owned-process guard checks host,
UID, PID start time and exact program/campaign arguments before every signal.
The [launch](launch.json) accounts for the geometry worker's remaining bound;
the joint future bound was 72.8% of available host RAM at launch.

The [resource-point driver](../../../../../scripts/validation/qualify_host_workspace.py)
must finish both gas/SMD hybrid energies and the geometry-level maximum-gradient
comparison against the preserved reference. Independent review is required
before adopting the execution profile for a completed geometry. This comparison
does not certify a full force vector, Hessian or conformer ensemble.
The independent CPU reference point on lambda0 continues separately.
Raw records and SCF checkpoints are in
`calcs/validation/2026-10-23_crest_refinement/baseline__design_pmdi_ethyl__occ_2_3_host_workspace_32GB/`.
