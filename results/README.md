# Results index

Use these stable paths. Completed quantities carry field-specific QC; blank values
are pending or excluded. Raw calculations and earlier source snapshots are retained.

| Result | Main entry | Current data |
| --- | --- | --- |
| Grafting free energies for all 75 candidates | [Visual report and structures](grafting/index.html) | [Method/records](grafting/README.md) · [Status](grafting/status.json) · [Table](grafting/current.csv) |
| Polymer-fragment and interaction models | [3D interaction report](polymer_models/index.html) | [Method/records](polymer_models/README.md) · [Interaction status](polymer_models/interaction_status.json) · [Graph inventory](polymer_models/current.csv) |
| Graft spacing and conditional grafting framework | [Exact occupancy cases and report](polymer_models/campaigns/2026-10-21_spacing_continuation/index.html) | [Reproduction guide](../docs/POLYMER_SPACING.md) · [Current pointer](polymer_models/spacing_current.json) · [Calculation status](polymer_models/campaigns/2026-10-21_spacing_continuation/status.json) |
| Candidate redox properties and comparisons | [Candidates](candidates/README.md) | [Workbook](candidates/latest.xlsx) · [All candidates](candidates/figures/data/all_candidates.csv) |
| Maintained figures and molecular structures | [Figure guide](candidates/figures/README.md) | [Comparison browser](candidates/figures/index.html) |
| Measured/calculated redox validation | [Validation guide](validation/README.md) | [Figure](figures/validation/matched_parent_validation.pdf) · [Actual values](figures/validation/matched_parent_validation.csv) |
| Candidate search, selection and source datasets | [Design/evaluation report](discovery/expansion/index.html) | [Discovery guide](discovery/README.md) · [Matched effects](discovery/expansion/matched_design_effects.pdf) · [Search inventory](discovery/systematic_families_2026-10-18/summary.json) |
| Full property-table arithmetic review | [Reviewed tables](validated/README.md) | [Dated review pointer](validated/current.json); later candidates use the active candidate feed |
| Research presentation | [Slides](../slides/README.md) | [Checked PDF](../slides/redox_workflow.pdf) |
| Detailed validation/provenance index | [Latest review](LATEST_REVIEW.json) | [Machine-readable result index](INDEX.json) |

The current synthesis route is grafting redox groups onto a preformed polymer.
See [the route definition](../docs/SYNTHESIS_ROUTE.md). The earlier vinyl-monomer
exploration is preserved separately under [vinyl_monomers](discovery/vinyl_monomers/README.md).

The [completed DBU/DIPEA pilot](grafting/base_assisted/index.html) supports retaining
both BQ designs and AQ-1-O as **base required**. Their modeled base-assisted
reactions are favorable; positive HCl-reference values alone do not exclude them.
Selectivity and yield remain untested. See the [current selection assessment](grafting/selection_assessment.json).

The [oligomer screen](polymer_models/README.md) targets weak intersite/segment
interaction. Strong persistent stacking and undesired intersite electron transfer
are design concerns. Its monomer/dimer/trimer pilot has no tetramers; neutral
electronic energies do not yet establish interaction free energies or transfer rates.
Open the [oligomer figure overview](polymer_models/summary.html) for construction,
checked 3D geometries and the completed electronic comparisons. The
[frozen-fragment decomposition](polymer_models/interaction_decomposition.json)
separates segment interaction from deformation as its new SPs qualify.
The [geometry audit](reviews/2026-10-20/polymer_geometry_recovery/geometry_audit.json)
found no flagged sanity issues in 27 starts and 17 completed DFT structures.
A separate [checkpointed ten-case recovery](polymer_models/campaigns/2026-10-20_geometry_recovery/plan.json)
addresses the original 100-step nonconvergences; use the current interaction status
for actual progress. All DFT geometry/energy work uses SMD acetonitrile.

The [configurable spacing framework](../docs/POLYMER_SPACING.md) now builds matched
adjacent/separated occupancy controls and keeps electronic nonadditivity separate
from conditional grafting ΔG. Its first phthalimide method-control panel has six
graphs and twelve qualified xTB warm starts, frozen for quantum evaluation.
Spacing DFT and native G are pending the resource-gated calculation stages;
prepared structures do not establish physical interaction strength or feasibility.

## Calculation records and provenance

Paths below are relative to the project root,
`/nfs/lambda_stor_01/homes/rzhu/0_redox` (also `/homes/rzhu/0_redox`).

| Record | Location |
| --- | --- |
| Candidate solution geometries, energies and thermal QC | `calcs/dft/<candidate_id>/<state>/` |
| Gas-phase UMA warm starts and descriptors | `calcs/uma/<candidate_id>/<state>/` |
| Exact active candidate snapshot and immutable export paths | [results/candidates/current.json](candidates/current.json) |
| Figure input snapshot, generated files and checksums | [results/candidates/figures/metadata/manifest.json](candidates/figures/metadata/manifest.json) |
| Experimental benchmark sources and conditions | [data/raw/validation/two_wave_mecn/README.md](../data/raw/validation/two_wave_mecn/README.md) |
| Scoped validation campaigns and their records | [results/validation/](validation/README.md) |
| Immutable table exports, checks and publication evidence | `results/reviews/`; follow [LATEST_REVIEW.json](LATEST_REVIEW.json) to the relevant scope |
| Grafting components and reaction identities | [results/grafting/species.json](grafting/species.json) · [reactions.json](grafting/reactions.json); [submitted campaign](grafting/campaigns/2026-10-19_mecn_reference/plan.json) with raw root `calcs/grafting/2026-10-19_mecn_reference/` |
| Oligomer graphs, occupancy and stereochemistry | [results/polymer_models/current.csv](polymer_models/current.csv); quantum calculations are pending |

The category status files own completion counts. A prepared structure is not a
computed result; an arithmetic check does not establish experimental accuracy.
Use the active candidate feed for current candidate coverage. The full-table
review is a dated cohort and can lag newly added structures.

## Updating results

Publish updates to the same current filenames, and keep earlier source snapshots
and evidence. Update the category guide/status and this index when a new result
type is added. Include the exact identities, method/solvent, pending quantities,
QC, raw paths and source manifest. Existing figure and slide refresh commands are
in their guides above. Root-level legacy CSV files are not the current entry point.
