# PMDI trimer force/energy diagnostic

The [diagnostic summary](summary.json) compares analytical solution gradients with
central energy differences at the saved, unqualified optimizer-returned geometry.
Both runs use r2SCAN-D4/def2-SVP/SMD acetonitrile and conservative SCF settings.
They make no geometry, energy or thermochemistry promotion.

At grid level 3, atom 7/Y (zero based) has analytic gradient +0.0001177801
Ha/Bohr while central differences give approximately −0.000063341 at both
0.0005 and 0.001 Bohr steps. This explains why more optimizer steps alone are
not an adequate repair. The level-5 run is internally consistent at its own
largest component, atom 14/Y. The original level-3 failing component still needs
the corresponding level-5 and implementation checks before claiming a correction.

The [frozen driver](driver_source.py), [launch plan](plan.json), resource preflight
and [small numerical evidence](evidence/grid3/status.json) preserve this finding.
The first grid-3 reservation was refused before computation; its successful
retry has a separate execution record. Raw displaced-coordinate SCFs remain in
`calcs/validation/polymer_gradient_diagnostic_2026-10-21/` with hashes in the summary.
