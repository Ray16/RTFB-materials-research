# Exact molar energy conversion

The maintained `redox.core.common.EV_KJ` now derives the eV-to-kJ/mol conversion
from the exact SI defining values of elementary charge and the Avogadro constant:
`1.602176634e-19 × 6.02214076e23 / 1000 = 96.48533212331002`.
`LEGACY_EV_KJ = 96.485` explicitly identifies historical-record arithmetic.
Conditional polymer reaction outputs and shared-reagent snapshots now report the
conversion factor alongside their values. Other maintained users of `EV_KJ`,
including reorganization and stability reports, inherit the corrected conversion.

[Independent decimal/cached-record checks](summary.json) and
[68 passing focused tests](tests.xml) cover exact-unit consistency, the cached
PMDI HCl/DBU reaction sums, balanced occupancy cycles, conformer partitions,
capacity and solution-reorganization regressions. The DBU cached reaction changes
by approximately −0.000355 kJ/mol relative to the historical conversion; its raw
eV components do not change. The independently calculated redox potential remains
in volts and is unaffected by this molar conversion.

No archived source, raw state or immutable manuscript snapshot is rewritten.
Historical snapshot checks continue to use their declared source factor; the
paper's primary reported values remain frozen. Namespace-isolated scientific
workers retain their original archived code. New maintained calculations use the
exact factor. This closes reviewer item m6 within that explicit scope; it does
not improve or validate electronic-structure accuracy.
