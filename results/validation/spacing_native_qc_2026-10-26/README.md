# Stored physical checks for the twelve spacing-control conformers

The [audit](summary.json) reconstructs each thermal receipt and checks every
hash-listed output before reading the exact native record and projected
partition. All twelve pose identities are retained. The largest Cartesian
force component ranges from 6.15046e-6 to 7.61831e-5 hartree/bohr, below the
unchanged 4.5e-4 gate. The table also reports Hessian asymmetry, the largest
sampled-column relative error and the lowest signed projected internal mode.

One pose has a −2.72392 cm^-1 mode, within the already declared −5..0 cm^-1
numerical tolerance. It enters qRRHO through its absolute value under the stated
policy; the table does not imply that every internal eigenvalue is positive.
No threshold was relaxed. These are stored checks of the original finite set,
not new force/Hessian evaluations or evidence of sampling convergence.

The manuscript and SI now include the actual values rather than only the
acceptance thresholds. The frozen artifact preserves all twelve rows and their
raw source paths/hashes. `scripts/validation/audit_spacing_native_qc.py` regenerates
the audit from the exact finite-state declaration in a new destination. The
source native records and their existing partitions were not modified.
