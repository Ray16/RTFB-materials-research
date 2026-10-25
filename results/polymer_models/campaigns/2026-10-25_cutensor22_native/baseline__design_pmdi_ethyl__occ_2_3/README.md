# Paired PMDI native thermochemistry: current controlled attempt

[Plan](plan.json) · [Actual startup](startup_verification.json) · [Live scientific stage](streamed_native/status.json) · [Runtime smoke](tensor_runtime/status.json) · [Controls](../../../../validation/cutensor22_full_spin_2026-10-25/comparison.json)

The exact 139-atom, two-graft PMDI model reuses its independently stationary geometry and canonical gas/SMD parent. The new attempt uses bounded density-fitting derivatives and an isolated CuPy 13.4.1/cuTENSOR 2.2.0 runtime, qualified separately against full singlet/doublet controls. Physical equations, precision and minimum/force gates are retained.

Pure-solute partial Hessian and H1 checkpoints have been saved under `calcs/polymer_models/2026-10-25_cutensor22_native/baseline__design_pmdi_ethyl__occ_2_3/streamed_execution/`. The subsequent orbital-response solve stopped with `Krylov subspace iterations diverge`. They remain **unqualified components**. Full response, projected Hessian/modes, minimum checks, native partition and independent reconstruction are required before publishing G. No paired G or neighboring-graft free energy is currently available.

All earlier allocation failures and the unsupported cuTENSOR 2.7 V100 smoke test remain at their original campaign paths. No failed geometry or partial tensor is promoted to a free-energy result.
