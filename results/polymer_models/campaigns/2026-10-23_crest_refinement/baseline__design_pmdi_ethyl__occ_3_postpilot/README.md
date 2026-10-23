# Exact site-3 baseline: continue after the checked resource pilot

The original bounded-workspace pilot passed at 76,862 MiB peak device memory.
Its full worker then stopped before optimization because the resource pilot had
already created the geometry-force evaluation directory. This branch uses fresh
full-calculation storage with the same task, geometry, method, optimizer and QC.
The original resource pilot, directory and failure receipt remain unchanged.

[plan.json](plan.json) binds the source and input hashes.
[parent_ownership_verification.json](parent_ownership_verification.json) checks
the original owner identities and absence of matching predictors before the
continuation is prepared. The [reused pilot receipt](energy/pilot/status.json)
retains its original process/timestamps and explicitly states that it is reused
evidence, not a new pilot calculation.

The normal H100 energy worker has started; its
[worker record](energy/workers/lambda13_gpu1_queue.json) identifies the actual
process and CPU assignment. Native G and the complete occupancy cycle remain
pending. The candidate observer now follows this branch and skips held task
claims before attempting additional launches.

Preparation: [prepare_polymer_postpilot_continuation.py](../../../../../scripts/pipeline/prepare_polymer_postpilot_continuation.py).
Raw calculations and durable geometry evaluations:
`calcs/polymer_models/2026-10-23_crest_refinement/baseline__design_pmdi_ethyl__occ_3_postpilot/`.
This remains a single-start feasibility branch; broader CREST convergence and
synthesis admission are separate requirements.
