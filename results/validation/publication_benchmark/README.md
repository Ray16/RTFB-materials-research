# Frozen retrospective redox benchmark

[Panel v1](panel_v1.json) freezes 39 strict experimental potential endpoints from
30 exact free-parent molecules before the remaining calculations complete.
All previously sourced strict endpoints are included. Lower-quality and excluded
source values remain in its source catalog. This is not a held-out test set.

The [current analysis](current.json) reports unrounded residuals, event- and
molecule-weighted errors, family/wave coverage, signed wave spacings, outliers,
and molecule/study cluster resampling. The current independent review admits 39/39 endpoints from 30 molecules and seven
original studies: MAE 0.247 V, RMSE 0.279 V and signed bias +0.190 V. All fixed endpoints are admitted. The last PMDI electronic-spin comparison passes with its alternative triplet 0.707 eV above the selected singlet; the alternative thermal spin gap remains uncomputed. Both molecule and study cluster intervals are
reported as descriptive resampling of this retrospective panel. No model
correction or outlier removal is fitted.

An additional [NDI experimental challenge](../ndi_external_challenge_2026-10-26/README.md)
has been identified separately. Its quasi-reversible peak midpoints are proxies,
its calculated states are pending, and it has not been added to this fixed panel.

The [latest independent review](../../reviews/2026-10-22/publication_benchmark_complete/summary.json)
and [reproduction guide](../../../docs/PUBLICATION_READINESS.md) distinguish
execution correctness, source quality, physical state QC and prediction accuracy.
The [native completion plan](native_completion_2026-10-21/plan.json) contains
48 initially missing native states in 22 molecule jobs. The separately held
BQ cache was independently found to be complete but unselected; its normal
thermal-policy selection was restored without recomputing a Hessian. Its [start verification](native_completion_2026-10-21/start_verification.json)
records actual processes, GPU residency, CPU affinity and thread caps.

Raw quantum records stay under `calcs/dft/<parent>/<state>/`. Native state
completion does not imply that a spin model, exact graft or polymer prediction
has been experimentally validated. Re-run the independent reviewer and checked
consumer before publishing a later snapshot.

The maintained [expanded benchmark figure](../../figures/validation/retrospective_redox_benchmark.pdf)
and [actual plotted values](../../figures/validation/retrospective_redox_benchmark.csv)
use the same independently checked, unrounded values. Page 7 of the updated
[slides](../../../slides/redox_workflow.pdf) presents this expanded panel; earlier
parent subsets retain their dated scope in Backup.
