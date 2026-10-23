# Bounded SMD execution: controls passed, large-case resource checks pending

The original 139-atom compact-build pilot later failed in `D*A` during solvent
response. Two 115-atom pilots completed point evaluations but exceeded their
80 GiB resource limit. Failures and their completed point evaluations remain
preserved; none establishes an optimized geometry or native G.

The successor evaluates `D diag(A) v` as `D (A v)` and its transpose as
`A (D^T v)`, avoiding a dense temporary. GPU4PySCF workspace estimates also
respect an 80 GiB allocation budget. Surface grids, basis, functionals, precision
and physical acceptance thresholds remain unchanged. The installed library is
never modified. The earlier qualified build-lifetime adapter is reused.

Both [singlet](controls/1/status.json) and [doublet](controls/2/status.json)
controls pass full original-versus-bounded energies, optimized coordinates,
gradients and Hessians. Maximum energy changes are 2.73e-12/9.09e-13 eV; largest
Hessian changes are 6.67e-12/8.98e-13 atomic units. Both exercised the new memory
queries and response contraction. This establishes numerical equivalence for
these controls; larger-system memory suitability remains a separate pilot.

[Launch/preflight](launch.json), [actual starts](actual_start_verification.json)
and [failed-owner death/capacity audit](failed_predictors_death_and_capacity.json)
record reservations, two disjoint CPU pairs and all native caps=2. The current
[PMDI recovery registry](../pmdi_occupancy_2026-10-23/baseline_recoveries.json)
selects three fresh same-model branches. Available preserved larger-case point
energies/gradients must agree before a new pilot is admitted. The current resource snapshot shows all eight suitable
H100s occupied; submission waits on the shared gate.
No smaller memory class or relaxed budget is substituted.
