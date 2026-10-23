# Solvent workspace and checkpoint qualification on 32 GB GPUs

[summary.json](summary.json) separates numerical control equivalence from large
model resource qualification. The calculations retain the existing solvent
surface, precision, integration grids and electronic-structure model.

The workspace adapter prebuilds persistent SMD storage before choosing density-
fitting workspaces, uses the earlier bounded contractions and saves each SCF
iteration plus a final HDF checkpoint on shared storage. Singlet and doublet
controls reproduce the original energies, coordinates, gradients and Hessians
within the existing tolerances. A separate bounded-gradient adapter also passes
both controls; its energy solver remains direct LU.

The isolated 115-atom PMDI starting-geometry pilot on lambda2 completed both
MeCN and gas-phase single-point SCFs, preserving their checkpoints. The later
geometry-level solvent gradient ran out of device memory while allocating its
dense derivative tensor. The failed attempt is preserved. Neither an optimized
state nor a native Hessian/free energy is qualified by these partial outputs.
Small-control success does not qualify the large-molecule gradient recovery.

The screened full- and range-separated density-fitting tensors require about
160.7 GB combined host storage; the conservative preflight required about
218.6 GB. Host RAM and disjoint CPU capacity therefore constrain suitable
workers in addition to device memory. The successful SCF checkpoints remain
under `calcs/validation/smd_workspace_32GB_2026-10-23/pilot_occ_2`.

Producers: [workspace adapter](../../../scripts/pipeline/smd_workspace_execution.py),
[fixed-geometry resource pilot](../../../scripts/validation/qualify_polymer_workspace.py),
[bounded-gradient adapter](../../../scripts/pipeline/smd_gradient_workspace.py).
The control receipts bind their exact source and installed-library hashes.
The large-gradient recovery and full optimization/thermal qualification remain open.
