# Fresh gas-phase reorganization reference replays

Two new GPU calculations exercise the corrected reference solver on the
B3LYP/6-31G(d,p), RI-J, neutral-singlet/cation-doublet four-point hole model.
Original saved records remain unchanged. This is an actual quantum recomputation,
not a retrospective assertion that old convergence metadata existed.

| Molecule | Saved legacy value | Fresh calculation | Curated computational reference |
| --- | ---: | ---: | ---: |
| Tetracene | 0.1114 | 0.1114 | 0.114 |
| Pentacene | 0.0931 | 0.0932 | 0.095 |

All values are eV. Changes are compared at the saved four-decimal precision.
The unrounded fresh values are 0.1114277458 and 0.0931529584 eV. Reference values
come from the existing [curated source table](../../reorg_anchors/oligoacene_inner.csv);
they are computational literature targets, not experimental measurements.
The two residuals do not establish general method accuracy or solution kinetics.

The [plan](plan.json), [preflight](lambda13_preflight.json), [launch](launch.json),
[actual starts](start_verification.json) and [terminal process check](terminal_verification.json)
record the system reservations, disjoint CPU pairs and actual worker ownership.
Both calculations finished with explicit SCF/optimizer convergence and full grid
response. The [independent review](review.json) passes 40 checks, reconstructing
the four-point sum and each directional term from four raw Hartree energies,
checking units/rounding and both optimized molecular graphs, and verifying all
legacy/source byte hashes. [Fresh records](calc/) retain energies and geometries.

The stated spin models and optimized geometries were reproduced; this work did
not add Hessian minimum proofs, a broader spin search or an experimental rate
comparison. The other legacy records in the 482-record inventory are not
retroactively certified. No production candidate value or canonical reference
record was replaced.
