# 64 GiB execution profile on H100: bounded gradients

Both full original-versus-modified implementation controls pass:
[singlet](controls/1/status.json) and [doublet](controls/2/status.json).
Maximum electronic-energy differences are below 6e-12 eV. Coordinates, gradients
and Hessians also meet the unchanged tolerances. The direct LU energy solver,
physical model, precision and solvent grids are retained; allocation/workspace
selection and bounded gradient contractions form the execution profile.

The [139-atom point campaign](pilot_occ_2_3/plan.json) follows the preserved
out-of-memory attempt. Its [launch](pilot_occ_2_3/launch.json),
[preflight](pilot_occ_2_3/preflight.json) and
[live resource status](pilot_occ_2_3/footprint_qualification.json) distinguish a
running qualification from a completed state. Screened density-fitting storage
requires a conservative 303.7 GB host-RAM budget, which fits the recorded H100
node preflight. The actual resource pilot must finish before suitability is
accepted. No optimized geometry, Hessian or free energy is claimed here.

[startup_verification.json](startup_verification.json) records the actual live
resource process and the parallel site-3 energy worker: separate GPUs, disjoint
CPU pairs and all native thread caps equal to two. The earlier V100 site-2
optimization continues separately. Existing scientific workers are not stopped.

Producer: [qualify_polymer_workspace_gradient.py](../../../scripts/validation/qualify_polymer_workspace_gradient.py),
which reuses the original tensor inventory, energy/gradient and continuity
checks with the separately controlled bounded-gradient adapter. All original
139-atom failed outputs remain at their source locations. New raw outputs:
`calcs/validation/smd_gradient_workspace_64GB_2026-10-23/`.
