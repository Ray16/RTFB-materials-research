# Residual-checked orbital-response execution

[Full spin and large-response review](comparison.json) · [CPU solver tests](../../../tests/test_residual_response_solver.py) · [Current full paired case](../../polymer_models/campaigns/2026-10-25_residual_native/baseline__design_pmdi_ethyl__occ_2_3/README.md)

The new solver leaves the supplied response operator and orbital-energy denominators unchanged. It solves the original matrix-free I+A equation with host-vector GMRES and GPU operator evaluations, then explicitly tests the true residual at a target 100 times tighter than the supplied absolute tolerance. No level shift is added. Bounded iteration limits and nonfinite/shape checks fail closed.

Full singlet/doublet control energies, stationary coordinates, gradients and Hessians are independently compared. Maximum Hessian differences are 9.09e-8 and 5.99e-8 atomic units, within the original 1e-6 bound. The large fixed-state probe restores exact saved MOs and pure-solute H1, adds the first-atom PCM nuclear derivative, and enables equilibrium solvent response just as in the full driver. Its independent unpreconditioned equation residual is 1.68e-10 hartree, below 1e-8. This is one atom's response, not a full molecular Hessian or free energy.

The first two probe setups are preserved: one allocated DF storage before building the PCM surface and exceeded memory; the next lacked the Hessian response object. The final solvated probe corrects the setup explicitly and does not repeat SCF. All scripts and source-bound plans remain distinct. Raw response arrays live under `calcs/validation/response_solver_2026-10-25/`; the validated runtime uses the existing isolated cuTENSOR 2.2 overlay.
