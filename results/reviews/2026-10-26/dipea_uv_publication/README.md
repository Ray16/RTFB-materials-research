# DIPEA UV and NMR source revision

[UV comparison](../../../validation/dipea_uv_anchor_2026-10-26/README.md) ·
[Independent source review](../../../validation/dipea_uv_anchor_2026-10-26/independent_review.json) ·
[Manuscript](../../../../manuscript/manuscript.pdf) ·
[Supporting Information](../../../../manuscript/supporting_information.pdf) ·
[Frozen evidence](../../../../manuscript/evidence/2026-10-26_dipea_uv/manifest.json) ·
[Slide 12](../../../../slides/colleague_workflow.pdf)

The checked primary sources report DIPEA log Kᴴ = 18.64(1) in MeCN at
25.0 ± 0.1 °C from UV titration against bromocresol green, and a separate
NMR estimate of 18.7 for DIPEAH hexafluorophosphate in MeCN-d3. Their
conditions and calibrations remain distinct; neither directly measures the
modeled chloride reaction. The manuscript uses the 25 °C source for the
conditional pKa-derived DIPEA/DBU and DIPEA/HCl comparisons, retains the NMR
estimate as an independent context, and explicitly avoids treating a hybrid
HCl-capture substitution as a validated full grafting result.

The current epoch extends the previous frozen NMR and numerical evidence;
older epochs and their figure files remain unchanged. The earlier publication's
55 source identities are independently [reconciled to retained bytes](historical_relocation.json)
where current presentation/manuscript paths have moved on. Focused tests reject
wrong solvents, values, source mismatches and promotion of a hybrid sensitivity
to a prediction. Manuscript and slide PDFs have separate compiler and visual
reviews, and the self-contained upload ZIP reproduces the checked text. No
Overleaf project has been connected or created.

The [full CPU regression](../../../validation/dipea_uv_anchor_2026-10-26/full_regression_summary.json)
passed 1,357 tests. Its 14 raw-package replay skips and one figure-environment
skip were run separately and passed. The later counterpoise coordinator tests
were also run separately; this software evidence does not certify physical
accuracy or a completed polymer ensemble.

Both third-party supporting-information PDFs and article XMLs are held in
`calcs/literature/2026-10-26_dipea_experiment/`, with source hashes; their
licenses are not replaced by the project code license. No new quantum jobs,
candidate admissions, grafting yields or ensemble-convergence claims follow
from this literature revision.
