# Checked occupied-state thermochemistry and resource pilot

The prospective 40-GiB host-memory pilot for `all_same__occ_4__conf003` has
completed both numerical/physical and resource reviews. The source-bound
collector verifies 2,449 bindings. This is one checked occupied conformation in
the alternate smooth-surface MeCN model; occupied ensembles and balanced
converged occupancy cycles remain incomplete.

Its free energy at the primary 50 cm−1 rotor cutoff is
−84296.34187317279 eV; the minimum projected frequency is +15.52 cm−1 and the
largest sampled Hessian-column discrepancy is 3.47e−6. The new result remains
separate from the paper's existing 75-state snapshot and from primary-model
PMDI/phthalimide cycles. Raw records and individual Hessian/SCF checkpoints are
preserved on shared storage.

The original 80-GiB inner guard and numerical equations were retained. A
prospectively declared outer guard imposed the stricter 40-GiB bound and retained
the 32-GiB free-memory floor. Tree-RSS peak was 25,618,284,544 bytes (23.86 GiB),
and process VmHWM was 25,600,720,896 bytes. GPU peak was 24,206 MiB. The tested
run passes; other conformers require an explicit matched constitution/state/basis
and their own prospective profile. The complete fleet is not retroactively
reclassified, and occupied GPUs remain unavailable.

The first collector stopped after the completed source-host review was not
visible at its declared shared path within 30 seconds. Exact owner checks later
confirmed that both scientific worker and collector had exited; the numerical
review file was present. A new source-host collector read the existing artifacts
and performed the separate resource review. No quantum calculation was relaunched,
no review was duplicated, and old ownership/status files remain unchanged.

[observer_status.json](observer_status.json),
[resource_review.json](resource_review.json) and
[native_review.json](native_review.json) preserve the completed outcome.
[profile.json](profile.json) records the prospective resource class.
The changing main dispatcher remains at
`results/validation/native_conveyor_2026-10-30/dispatch/status.json`; successful
electronic evaluation is never substituted for missing thermochemistry.
