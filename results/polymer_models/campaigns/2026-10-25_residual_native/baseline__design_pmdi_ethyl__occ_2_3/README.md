# Exact paired PMDI: full native attempt with residual-checked response

[Plan](plan.json) · [Actual startup](startup_verification.json) · [Live stage](streamed_native/status.json) · [Controls and first-atom diagnostic](../../../../validation/response_solver_2026-10-25/comparison.json)

The new attempt preserves the exact stationary geometry, electronic parent, grids, functional/basis and all native physical checks. Its numerical response solver uses normalized GMRES with an explicit original-equation residual check; no level shift or looser tolerance is introduced. Complete singlet/doublet Hessians reproduce the original path within its predeclared tolerance, and the actual 139-atom first-atom solvent response passes an independent unpreconditioned residual test.

This is **a running full-Hessian attempt**, not a completed free energy. The ordinary projected modes, minimum, finite-difference, nuclear-partition and independent reconstruction checks remain required. Raw partial derivatives and SCF checkpoints persist under `calcs/polymer_models/2026-10-25_residual_native/`. The earlier allocation and Krylov failures remain at their original campaign paths.
