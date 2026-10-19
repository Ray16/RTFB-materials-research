# Rigid-motion exclusion in native thermochemistry

The [mode-count audit](mode_count_audit.json) found extra modes in all 232
inspected current-geometry native records. Their stored PySCF spectra exclude
translation and rotation, but their xTB thermal logs did not. xTB 6.7.1's
[Turbomole input path](https://github.com/grimme-lab/xtb/blob/v6.7.1/src/prog/thermo.f90#L102-L107)
mass-weights the input and diagonalizes it without the rigid projection that its
DFTB+/ORCA paths perform. Previous arithmetic audits reproduced these logs and
therefore did not prove physically correct thermal partitioning.

`src/redox/qm/vibrations.py:project_rigid_modes` now provides explicit
mass-weighted projection with SVD, preserving internal negative curvature.
Eight regression cases check rigid nullspaces, linear/nonlinear molecules,
idempotence, coordinate translation/scaling and internal stretching/instability.
This helper is not yet integrated into the live native producer or primary readers.

The [bounded CPU replay](projected_replay/summary.json) reuses frozen Cartesian
Hessians, writes projected Cartesian xTB inputs, and recomputes 25/50/100 cm⁻¹
qRRHO corrections. It preserves canonical data and completed replay chunks.
An [independent review](independent_review/summary.json) reconstructs the projector
using reorthogonalized Gram-Schmidt, checks serialization and all mode counts,
replays thermal arithmetic and spin partitioning, and verifies input hashes.
It passes **7,854 checks for 231 records**, with zero failed review checks.

One PMDI methoxyethoxyethyl neutral record remains excluded: its genuine
0.391 cm⁻¹ internal mode is below xTB's retained-frequency threshold. No threshold
was relaxed or missing mode silently accepted.

The [candidate impact diagnostic](candidate_potential_impact.csv) and
[summary](candidate_impact_summary.json) recompute the matched Fc effect for the
16 previously complete additional ladders. Their largest individual potential
shift is 0.0493 V. These are diagnostic results, not newly admitted predictions.
Fc itself shifts by +0.007984 V; its projected diagnostic absolute value is
4.379196 V. Controls and spin-ensemble decisions still need consistent treatment.

The [publication hold](../../../config/thermal_publication_review.json) keeps
affected current candidate fields blank. Historical snapshots and slides preserve
their earlier numbers with explicit review labels. Release requires producer and
both-reader integration, independently verified matched references and exports,
and re-evaluation of the shortlist and experimental benchmark. Live numerical
producer files remain frozen while existing source-bound calculations run.
