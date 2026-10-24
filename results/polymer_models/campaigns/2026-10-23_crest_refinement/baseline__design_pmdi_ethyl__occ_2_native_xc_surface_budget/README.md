# Site-2 native thermochemistry after explicit force-validation refinement

The energy parent is the [separate four-step qualification capsule](../baseline__design_pmdi_ethyl__occ_2_fine_force_qualified/README.md), which preserves the original failed coarse check and unchanged force-error tolerance. [plan.json](plan.json) and the [energy receipt](energy/done/baseline__design_pmdi_ethyl__occ_2__pose0.json) retain its exact original geometry/active energy and named validation procedure. Preparation creates no G and repeats no optimization.

The native resource pilot has an actual reserved V100 start on lambda2:
[verified process/GPU/CPU resources](../../../../validation/pmdi_site2_fine_force_2026-10-23/native_startup_verification.json). Both CPU cores are disjoint from the observer; all six native caps are two. The measured-cavity host forecast is 214.93 GiB, retaining the original cache and headroom formula.

No accepted native G exists at this checkpoint. The original native DFT Hessian, symmetry, finite-gradient, stationary-gradient, imaginary-mode, spin and partition checks remain required. This pilot runs in separate raw storage; it cannot alter the old failed source or manufacture a conformer ensemble. Raw paths, controls and source hashes are given in the plan. Exact-candidate polymer cycles and all-axis eligibility remain pending.
