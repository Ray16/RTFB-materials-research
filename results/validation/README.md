# Validation results

Return to the [results index](../README.md). Validation is scoped by property and
cohort; the available experimental benchmark does not validate every property or
the behavior of grafted polymers.

| Question | Authoritative record |
| --- | --- |
| How do absolute calculated potentials compare with measurements? | [Matched-parent benchmark](../reviews/2026-10-19/full_matched_parent_final/README.md): 12 potentials from eight free parents; [measured/calculated values](../figures/validation/matched_parent_validation.csv) and [comparison figure](../figures/validation/matched_parent_validation.pdf) |
| What are the observed errors? | [Benchmark summary](../reviews/2026-10-19/full_matched_parent_final/summary.json): MAE 0.209 V, RMSE 0.246 V; the phenothiazine second-wave outlier remains included |
| Is the corrected thermal partition implemented and checked? | [Thermochemistry qualification](rigid_mode_thermal_2026-10-18/README.md) and [publication admission](../../config/thermal_publication_review.json) |
| Are assembled property tables arithmetically consistent? | [Dated full-table review](../validated/README.md); independent consistency checks are separate from predictive accuracy |
| Which properties can enter candidate selection? | [Model limits](../../config/model_validation_limits.json) and the field-specific QC in the [active candidate feed](../candidates/current.json) |
| Where are the latest scoped checks and remaining limits? | [Latest review pointer](../LATEST_REVIEW.json) |

The benchmark's exact sources and experimental conditions are recorded in
[comparisons.csv](../reviews/2026-10-19/full_matched_parent_final/comparisons.csv)
and the [experimental data guide](../../data/raw/validation/two_wave_mecn/README.md).
No fitted offset is transferred to a grafted candidate. Exact-graft experimental
accuracy, solution kinetics and polymer performance remain separate open tasks.

Individual campaign directories retain plans, ownership/launch records, logs,
raw replay artifacts and independent reviews. Their dated counts describe that
campaign. Current candidate calculations live in `calcs/dft/<id>/<state>/` from
the project root; follow a review's manifest for isolated alternatives and repairs.
