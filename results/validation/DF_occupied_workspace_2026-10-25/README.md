# Bounded density-fitting derivatives and tensor-runtime controls

[Two-spin/cap controls](caps/) · [Compatible cuTENSOR comparison](../cutensor22_full_spin_2026-10-25/comparison.json) · [Runtime guide](../../../environments/cutensor22/README.md)

The versioned adapters batch auxiliary and occupied-orbital work while retaining all atoms, orbitals and scientific settings. Complete singlet/doublet Hessians are compared against the original path at two occupied-block caps. Partial Hessian and H1 arrays are durable **unqualified** checkpoints, bound to coordinates, molecular orbitals, basis and SCF energy. Completing a partial derivative is not completing a Hessian or free energy.

The large paired PMDI case exceeded memory at successive original allocation stages. Those attempts remain preserved. The current V100-compatible runtime uses the separately controlled cuTENSOR 2.2 overlay and has checkpointed both pure-solute partial Hessian and H1. The subsequent Krylov orbital-response solve diverged; no G was published. Full response, projected modes and physical minimum/partition checks remain required. See the [current case](../../polymer_models/campaigns/2026-10-25_cutensor22_native/baseline__design_pmdi_ethyl__occ_2_3/streamed_native/status.json).

Raw large tensors stay under `calcs/polymer_models/2026-10-25_cutensor22_native/`; isolated fixed-orbital comparisons under `calcs/validation/DF_occupied_large_2026-10-25/` and `calcs/validation/DF_cutensor22_large_2026-10-25/` are subroutine controls only. Full-molecule 16 GB execution is not admitted.
