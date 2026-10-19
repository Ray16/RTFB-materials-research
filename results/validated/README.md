# Reviewed property tables

**Corrected property-table arithmetic review, 2026-10-19.** These tables use the
independently qualified internal-mode thermochemistry and matched Fc reference.
The earlier unprojected review remains preserved under
`results/reviews/2026-10-18/post_finished_spin_reports/` and is superseded.

The [current review](current/) contains the full property-table assembly and
its independent arithmetic audit. It passed **32,492 checks, zero discrepancies**,
with unchanged inputs during the review and a second inventory check before
publication. These are dated, immutable tables; later
calculations may supersede individual fields.

- [Redox potentials](current/redox_potentials.csv)
- [Scorecard](current/scorecard.csv)
- [Computed workbook](current/Candidates_computed.xlsx)
- [Review status](current/review_status.json) and [audit results](current/validation/audit_recompute.csv)

Arithmetic consistency is separate from physical and experimental validation.
Pending or excluded fields and spin confidence retain their own QC. This review
does not establish experimental accuracy or validate solution kinetics.

For the actively computed additional candidates, use the more frequently refreshed
[candidate results](../candidates/README.md). Earlier root-level CSVs and prior
review snapshots remain historical evidence.

The current review also applies the evidence-bound phenothiazine second-oxidation
selection limit. Calculated values remain recorded; affected window/capacity
objectives cannot certify primary shortlist membership. See
[`config/model_validation_limits.json`](../../config/model_validation_limits.json).

The current snapshot includes both new cata-imide ladders and gas-diagnostic QC.
The repaired cata-ADI neutral uses a reviewed finer CDS derivative step with the
original geometry and active energy unchanged. Per-record external qualification
evidence is now fingerprinted, so a changed gradient artifact cannot hide behind
an unchanged native-record JSON. The software suite passes 544 tests.

The current snapshot also exports computational-only lambda uncertainty evidence
and guards its prior floor. See the [review](../validation/uncertainty_scope_2026-10-19/README.md).
Current numerical property/scorecard/Pareto tables remain unchanged.
