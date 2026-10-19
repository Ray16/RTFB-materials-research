# Hexyl-NDI contact diagnosis and distinct seed preparation

The preserved radical geometries remain unqualified. This CPU-only follow-up
does not repeat their DFT/Hessians, alter native thresholds, change the installed
solver or replace canonical records.

## Contact diagnosis

The old derivative discrepancy involves H atoms 40 and 71 (zero-based). Their
separations are 3.200012 and 3.199992 A in the two frozen geometries. The tagged
upstream CDS implementation uses a 1.2-A hydrogen radius plus a 0.4-A solvent
radius, giving an expanded H/H sphere-contact distance near 3.2 A.[^1]

The first [installed-library sweep](contact_sweep/report.json) located a force
jump near that contact. A [precision-resolved sweep](contact_sweep_refined/report.json)
centers the contact at 3.2000000953674315 legacy A, consistent with the source's
unsuffixed radius literal. This is a source/behavior correspondence; the build
identity of the installed binary and tagged source is not claimed.

| Diagnostic | Minus geometry | Plus geometry |
| --- | ---: | ---: |
| Limiting full CDS gradient jump (Ha/Bohr) | 0.000219830 | 0.000219831 |
| Curvature growth per tenfold smaller cross-contact step | 10x | 10x |

The [independent review](contact_independent_review.json) passes 38 checks. The
jump divided by the original `2h` predicts a 0.219829879-Ha/Bohr2 residual,
compared with 0.219800296 from the preserved full quantum gradient/Hessian
column: a relative difference of 0.01346%. The contact therefore accounts for
the observed residual scale. Both failed geometries lie very close to this
non-smooth cavity contribution. Taking a smaller step on one side cannot certify
the original cross-contact derivative test. No Hessian is symmetrized or admitted.

## New starting conformers

All ten original FF conformers were preoptimized using the same GFN2-xTB/ALPB(MeCN)
charge -1, doublet model. The declared [screen](plan.json) checks exact H topology,
heavy-atom RMSD, changed chain torsions and host-only CDS step consistency.
Both preserved CDS control matrices reproduce exactly before screening.

The [independent preparation review](independent_review.json) passes 149 checks,
including separate alignment/torsion calculations, energy parsing, raw matrix
checks and the selection rule. The two selected starting conformers are:

| Seed | Heavy RMSD from canonical geometry (A) | Chain torsions changed by at least 60 degrees |
| --- | ---: | ---: |
| [conf_08](conf_08/xtbopt.xyz) | 4.307 | 6 |
| [conf_09](conf_09/xtbopt.xyz) | 5.029 | 7 |

These are ALPB-preoptimized inputs, not verified DFT minima. The lowest ALPB-energy
input did not pass the CDS screening heuristic; all ten records remain available.
Screening does not establish native energy ordering or an equilibrium partition.
A higher-energy conformer must not replace an unqualified lower minimum merely
because its derivatives pass. Future DFT work must compare energies and conformer
scope, reject collapse back into the failed basin, and independently qualify any
native partition before considering adoption. Canonical potentials stay pending.

Raw CDS matrices and gradient arrays remain on shared storage under this folder;
their producing inputs and source hashes are recorded in the JSON evidence.
The first review attempt used an incorrect ASE import location; its log is
preserved. Correcting that import repeated no molecular calculation.

## Private native-level continuation

The [predeclared native plan](native_seed_plan.json) now advances `conf_08` and
`conf_09` to separate strict r2SCAN-D4/def2-SVPD/SMD(MeCN) optimizations. Both
workers actually started on reserved V100-32GB GPUs on lambda2, with disjoint CPU
pairs and two native threads; see [start verification](native_start_verification.json).
Each seed has one attempt and a 120-step limit. The active wB97M-V/TZVPD energies
are computed before the standard native Hessian and corrected partition, retaining
electronic comparisons if derivative qualification subsequently fails.

All outputs live under `native_seeds/`; canonical source hashes are frozen and no
promotion is performed. Both final basins need independent hydrogen/integrity,
collapse, energy-ordering and native-partition review against each other and the
preserved lower unqualified basin. These calculation starts establish neither
native minima nor a completed potential pair. The physical thresholds and original
failed contact evidence remain intact.

[^1]: GPU4PySCF v1.8.1, [mnsol.F](https://github.com/pyscf/gpu4pyscf/blob/v1.8.1/gpu4pyscf/lib/solvent/mnsol.F), solvent-radius setup and `VDWRAD`. Source URL, blob identity and checksum: [upstream_source_metadata.json](upstream_source_metadata.json). The downloaded source is retained locally for inspection; it is not substituted into the installed library.
