# Signed segment deformation and actual relaxation

The three electronic decompositions reproduce from the six frozen-fragment
single points and their independently optimized monomer references. No new
calculation or geometry change was made by this audit.

| Pair | Pair optimization steps | Pair relaxation RMSD (Å) | Fragment contributions (kJ/mol) | Sum (kJ/mol) |
| --- | ---: | ---: | ---: | ---: |
| Phthalimide, parallel | 76 | 0.700 | −0.680, +1.241 | +0.561 |
| Phthalimide, antiparallel | 57 | 0.245 | +0.142, +0.351 | +0.493 |
| PMDI, antiparallel | 93 | 0.501 | +1.002, −0.620 | +0.382 |

Small sums therefore do not establish that the complexes barely relaxed. Two
sums contain cancellation. Individual frozen-to-relaxed monomer RMSDs span
0.066–0.612 Å. Each segment contains 35 or 45 atoms, including its ethylbenzyl
backbone/cap model; these are not long polymer chains or bare redox cores.

The deformation term is a **signed difference from one sampled monomer
reference**, evaluated with a higher-level single-point protocol than the geometry
optimizer. A globally minimizing monomer reference is not established, and the
term is not a positive-definite strain penalty. The stored total decomposition
and its operational definition remain valid; their interpretation needs this
qualification. The historical arithmetic used 96.485 kJ/mol/eV. The new report
also records exact-SI conversion; original values are preserved.

Saved optimizer/SCF convergence and independent stationary-gradient records are
checked. This does not reconstruct the full constraint history or establish
Hessian minima, ensemble convergence, BSSE corrections, pure-core binding or
electron-transfer rates. The [source-bound summary](summary.json) and
[producer snapshot](producer_source.py) record the checked evidence.
