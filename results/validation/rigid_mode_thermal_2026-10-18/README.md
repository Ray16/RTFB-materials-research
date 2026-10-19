# Rigid-motion exclusion in native thermochemistry

The corrected thermal partition is integrated and independently qualified.
Current candidate exports require it; raw native Hessians and historical logs
remain preserved. [The publication review](../../../config/thermal_publication_review.json)
releases the implementation hold while retaining per-field physical QC.

## Cause and correction

The [mode-count audit](mode_count_audit.json) found extra modes in all 232
inspected current-geometry records. PySCF's minimum checks excluded rigid motions,
but the Turbomole xTB thermochemistry path did not. Earlier arithmetic reviews
reproduced these logs without detecting the physical partition error.

`qm.qrrho` constructs the internal mass-weighted subspace and passes its explicit
frequencies to the installed xTB thermodynamic kernel. Genuine soft internal modes,
including the 0.391 cm⁻¹ PMDI mode, are retained. `qm.xtb_kernel` isolates that kernel
in a dedicated CPU process because importing it beside PySCF's D4 libraries caused
a symbol conflict. No installed library or acceptance threshold was changed.
`qm.projected_thermal` writes source-addressed corrected caches alongside raw data.
The native producer, primary reader and independently implemented audit require
these caches. `redox thermal-postprocess` reuses qualified existing Hessians on CPU.

## Authoritative evidence

- [Sealed qualification](../../../config/projected_thermal_validation.json):
  232 cases, **6,961 checks, zero failures**. The independent analytic equations
  and installed isolated kernel agree within 7.86e−12 eV. The qualification also
  preloads PySCF/D4 in the parent process to exercise the actual library conflict.
- [Current-cache review](isolated_current_independent_review/summary.json):
  236 records, **3,787 checks, zero failures**. This is a fixed reviewed cohort;
  later states require their own source-addressed partition and field QC.
- [Matched-parent review](../../reviews/2026-10-18/matched_parent_projected_checked/summary.json):
  seven potentials, four parent molecules, three studies, **171 checks, zero
  discrepancies**; RMSE 0.1714 V, MAE 0.1471 V. Exact-graft accuracy remains unvalidated.
- Corrected matched Fc absolute reference: **4.3791959016 V**. The older
  4.3712115520 V and unprojected benchmark statistics are superseded.
- [Current candidate feed](../../candidates/current.json) and
  [candidate workbook](../../candidates/latest.xlsx): field-specific release,
  blank pending values, preserved exact-structure identities and matched controls.

The [corrected global table assembly](../../reviews/2026-10-19/corrected_geometry_driver_tables/review_status.json)
passes 28,967 arithmetic checks with stable inputs and zero discrepancies.
Solution kinetics and exact-graft experimental validation remain incomplete.

## Preserved investigation history

The [first projected CLI replay](projected_replay/summary.json) and its
[independent review](independent_review/summary.json) passed 7,854 checks for
231 records but lost one real soft internal mode to the CLI's frequency filter.
That method was diagnostic and was superseded by the explicit-mode kernel.
The [earlier status](README_before_integration.md), source snapshots, original
logs and candidate impact diagnostic remain available for provenance.
