# Recovered final checks on the exact two-graft geometry

The preceding optimizer converged and saved its 139-atom coordinates. A later
displaced-energy force check exhausted GPU memory. This [recovery](plan.json)
reuses the exact optimizer output and repeats the original final-force checks
with collection before mean-field construction. The small singlet/doublet
controls preserve energies, gradients and full Hessians to numerical precision.

The [completed geometry report](geometry/status.json) passes: maximum gradient
7.62e-6 hartree/bohr, maximum energy-difference force error 1.03e-7 hartree/bohr
against the unchanged 1e-6 limit, and 27,894 MiB observed GPU peak. Graph,
stereochemistry and the intended open backbone also pass. No optimization was
repeated. All earlier evaluations and the interrupted force check remain intact.

This establishes a stationary geometry, not a Hessian minimum or active hybrid
energy. The latter proceeds in its separately checkpointed canonical SP stage.
Raw force records and SCF files are under
`calcs/polymer_models/2026-10-24_geometry_finish/`.
