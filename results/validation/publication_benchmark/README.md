# Frozen retrospective redox benchmark

[Panel v1](panel_v1.json) freezes 39 strict experimental potential endpoints from
30 exact free-parent molecules before the remaining calculations complete.
All previously sourced strict endpoints are included. Lower-quality and excluded
source values remain in its source catalog. This is not a held-out test set.

The [current analysis](current.json) reports unrounded residuals, event- and
molecule-weighted errors, family/wave coverage, signed wave spacings, outliers,
and molecule/study cluster resampling. Missing values remain pending. Four studies
currently contribute admitted values, so the declared five-cluster rule leaves
the study-level interval blank. No model correction or outlier removal is fitted.

The [independent source review](../../reviews/2026-10-21/publication_benchmark_coverage/summary.json)
and [reproduction guide](../../../docs/PUBLICATION_READINESS.md) distinguish
execution correctness, source quality, physical state QC and prediction accuracy.
The [native completion plan](native_completion_2026-10-21/plan.json) contains
48 missing native states in 22 molecule jobs and one explicitly held partial
attempt. Its [start verification](native_completion_2026-10-21/start_verification.json)
records actual processes, GPU residency, CPU affinity and thread caps.

Raw quantum records stay under `calcs/dft/<parent>/<state>/`. Native state
completion does not imply that a spin model, exact graft or polymer prediction
has been experimentally validated. Re-run the independent reviewer and checked
consumer before publishing a later snapshot.
