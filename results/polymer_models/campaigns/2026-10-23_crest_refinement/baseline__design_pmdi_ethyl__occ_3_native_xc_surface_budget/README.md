# Exact PMDI site-3 native thermochemistry pilot

The parent electronic workflow passed stationary geometry, force, active-energy
and connectivity checks. The [energy receipt](energy/done/baseline__design_pmdi_ethyl__occ_3__pose0.json)
points to that original report without recomputation or a fabricated thermal
correction. [plan.json](plan.json) preserves task, model, geometry and source hashes.

A normal native-thermal resource pilot was launched on a reserved V100 GPU on
lambda1 with two disjoint CPU cores and all native thread caps equal to two.
The checked execution controls are indexed in
[the control guide](../../../../validation/xc_workspace_lifetime_28GiB_2026-10-23/README.md).
Host admission uses the measured cavity-size forecast of 216.16 GiB, retaining
the original 128 GiB cache limit and extra 128 GiB headroom. The earlier
preparations that failed admission were never scientific starts and are preserved.

No accepted native G has been produced. The native DFT Hessian, sampled derivative
checks, physical mode checks and nuclear/spin partition reconstruction remain
required. A resource pilot or checked electronic energy cannot fill the missing
G field. Broader conformer sampling and candidate route/coupling/accessibility
requirements remain independent.

The V100 pilot subsequently failed in the analytic PCM Hessian response with a device-memory allocation error. Its checkpoints and failure records are preserved; no G was accepted. A [fresh H100 qualification](../baseline__design_pmdi_ethyl__occ_3_native_xc_H100/README.md) retains the same electronic parent and physical gates in separate storage.
