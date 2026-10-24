# Between-SCF memory collection

Both full singlet and doublet dense-reference controls pass after adding a
collection step before construction of each mean field. The largest
energy difference is 1.82e-12 eV and full-Hessian difference is 2.93e-12 atomic
units. Geometry and gradient checks pass the original bounds. Installed
libraries, physical equations, numerical grids and thresholds are unchanged.

The [adapter](../../../scripts/pipeline/smd_gc_workspace.py) collects unreachable
objects and frees unused GPU-pool blocks. Reachable scanners and densities stay
intact. The [qualifier](../../../scripts/validation/qualify_gc_workspace.py)
requires both exercised full-spin controls. Source hashes, small matrices and
launch/ownership records are retained in this directory.

The exact two-graft optimization returned convergence and saved its coordinates.
Its subsequent displaced-energy force check exhausted GPU memory. A separate
recovery reuses that optimizer result and repeats the original final force
checks with the lifecycle hook. Its source is
`scripts/pipeline/finish_polymer_geometry.py`; native G and the matched cycle
require subsequent stages. Passing these small controls alone does not qualify
the returned polymer geometry.
