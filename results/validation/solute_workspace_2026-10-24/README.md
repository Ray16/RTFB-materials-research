# Pure-solute workspace controls

Both full singlet/doublet controls pass original-versus-adapted optimized
coordinates, gas/SMD MeCN energies, gradients and complete Hessians. The largest
Hessian difference is 1.03e-11 atomic units, below the original 1e-6 limit.
[summary.json](summary.json) records the independent operator and residency
checks. The [software regression](regression.xml) has 1,125 passes and one skip;
that figure test also [passes in its separate environment](figure_test.xml).

The profile composes the checked streamed/weighted solvent response with
temporary host storage of SMD `S`, `D` and `K_LU` matrices during the pure-solute
`partial_hess_elec` and `make_h1` calls. Access to a suspended matrix raises an
error. All matrices return with identical FP64 bits before solvent derivatives
or response run. Both controls exercise two suspensions and six exact matrix
round trips. Installed libraries and prior source-bound producers are unchanged.

The [previous full native attempt](../../polymer_models/campaigns/2026-10-24_weighted_native/baseline__design_pmdi_ethyl__occ_2_3/failure_review.json)
failed in a solute density-fitting contraction before the weighted response was
called. Its SCF checkpoint is preserved. This does not invalidate the independently
checked geometry, canonical gas/SMD energies or the passed fixed-density
solvent-response test.

The [fresh native branch](../../polymer_models/campaigns/2026-10-24_solute_workspace_native/baseline__design_pmdi_ethyl__occ_2_3/plan.json)
retains the exact electronic parent and sampled response proof. Its
[status](../../polymer_models/campaigns/2026-10-24_solute_workspace_native/baseline__design_pmdi_ethyl__occ_2_3/streamed_native/status.json)
owns completion. Raw SCF/Hessian work is under
`calcs/polymer_models/2026-10-24_solute_workspace_native/`.
Large-state resource, minimum and free-energy qualification remain pending.
