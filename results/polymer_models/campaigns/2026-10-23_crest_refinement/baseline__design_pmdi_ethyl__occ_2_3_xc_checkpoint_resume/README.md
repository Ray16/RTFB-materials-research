# Exact PMDI two-graft restart from a checked evaluated checkpoint

[plan.json](plan.json) defines the same neutral-singlet 139-atom graph and physical
model as the preserved failed full optimizer. Evaluated frame 23 supplies only
starting coordinates. It is not relabeled as a stationary state or free energy.
Its graph, stereochemistry and sampled total-energy derivatives passed
[separate checkpoint checks](../../../../validation/xc_workspace_lifetime_2026-10-23/README.md).
The full optimizer now uses the controlled XC temporary-lifetime adapter.

The original same-graph hybrid-resource point and the new checkpoint-gradient
test are separate evidence: no new hybrid SCF is claimed for preparation, and
no prior energy is copied into a final-state slot. Full optimization, exact final
energy, physical geometry/force checks, normal native thermal qualification and
ensemble assessment remain required.

The candidate observer registered this branch for the H100 node. At the documented
checkpoint, submission observers are live but no scientific worker has started:
all node GPUs are occupied by resident processes. Resource contention does not
consume a molecular attempt. Launch handles and claims are checked before any
new start. The original failed branch and all raw evaluations are preserved.

Small result records live here; raw calculation storage is specified by the
plan. The source baseline remains a one-start control, independent from the
three broader PMDI CREST source searches. Neither baseline preparation nor
successful resource qualification admits this candidate for synthesis.

The optimizer subsequently started on an H100 and saved 23 complete evaluated checkpoints. A delayed foreign CUDA process was found on that device; only this owned worker was stopped. Its partial trajectory and peer process were preserved. The [exclusive continuation](../baseline__design_pmdi_ethyl__occ_2_3_exclusive_resume/README.md) uses the selected evaluated frame as an unqualified warm start. No source E or G was admitted from the interrupted run.
