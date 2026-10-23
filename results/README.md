# Results index

Use these stable paths. Completed quantities carry field-specific QC; blank values
are pending or excluded. Raw calculations and earlier source snapshots are retained.

| Result | Main entry | Current data |
| --- | --- | --- |
| Grafting free energies for all 75 candidates | [Visual report and structures](grafting/index.html) | [Method/records](grafting/README.md) · [Status](grafting/status.json) · [Table](grafting/current.csv) |
| Polymer-fragment and interaction models | [3D interaction report](polymer_models/index.html) | [Method/records](polymer_models/README.md) · [Interaction status](polymer_models/interaction_status.json) · [Graph inventory](polymer_models/current.csv) |
| Graft spacing and conditional grafting framework | [Exact occupancy cases and report](polymer_models/campaigns/2026-10-21_spacing_cds_recovery/index.html) | [Reproduction guide](../docs/POLYMER_SPACING.md) · [Current pointer](polymer_models/spacing_current.json) · [Independent validation](polymer_models/workflow_validation.json) |
| Candidate redox properties and comparisons | [Candidates](candidates/README.md) | [Workbook](candidates/latest.xlsx) · [All candidates](candidates/figures/data/all_candidates.csv) |
| Maintained figures and molecular structures | [Figure guide](candidates/figures/README.md) | [Comparison browser](candidates/figures/index.html) |
| Measured/calculated redox validation | [Validation guide](validation/README.md) | [Figure](figures/validation/matched_parent_validation.pdf) · [Actual values](figures/validation/matched_parent_validation.csv) |
| Publication preparation and benchmark completion | [Readiness guide](publication_readiness/README.md) | [Requirements](publication_readiness/current.json) · [Frozen benchmark](validation/publication_benchmark/README.md) |
| Candidate search, selection and source datasets | [Design/evaluation report](discovery/expansion/index.html) | [Discovery guide](discovery/README.md) · [Matched effects](discovery/expansion/matched_design_effects.pdf) · [Search inventory](discovery/systematic_families_2026-10-18/summary.json) |
| Full property-table arithmetic review | [Reviewed tables](validated/README.md) | [Dated review pointer](validated/current.json); later candidates use the active candidate feed |
| Working in silico manuscript | [Draft PDF](../manuscript/manuscript.pdf) | [Source and build](../manuscript/README.md) · [Early review questions](../manuscript/review_notes.md) · [Frozen evidence](../manuscript/evidence/2026-10-22_review_checked/manifest.json) |
| Research presentation | [Slides](../slides/README.md) | [Checked PDF](../slides/redox_workflow.pdf) |
| Detailed validation/provenance index | [Latest review](LATEST_REVIEW.json) | [Machine-readable result index](INDEX.json) |

The current synthesis route is grafting redox groups onto a preformed polymer.
See [the route definition](../docs/SYNTHESIS_ROUTE.md). The earlier vinyl-monomer
exploration is preserved separately under [vinyl_monomers](discovery/vinyl_monomers/README.md).

The [completed DBU/DIPEA pilot](grafting/base_assisted/index.html) supports retaining
both BQ designs and AQ-1-O as **base required**. These signs are conditional on the computed base term; the DBU capture check differs from the MeCN pKa scale, and DIPEA lacks an exact checked anchor. positive HCl-reference values alone do not exclude them.
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
The spacing panel has all 12 checked electronic conformers and all six selected
native free energies. The [conditional-grafting figure](figures/polymer_models/conditional_grafting.pdf)
shows neighboring-site ΔΔG of −5.15 kJ/mol (adjacent) and −11.23 kJ/mol (one spacer).
These include backbone, solvent and conformational response. Rotor-cutoff
sensitivity retains the signs but does not measure conformational uncertainty.
The [independent cycle review](polymer_models/workflow_validation.json) includes
DBU-assisted conditional grafting values and checks atom/charge balance and cycle
closure. The original oligomer panel is 26/27 checked structures; the final trimer
returned optimizer convergence but failed final force/energy derivative validation and remains excluded. Broader motif sensitivity, conformer sampling and full workflow closure remain open.
Follow [live completion status](polymer_models/completion/status.json).

The [retrospective benchmark](validation/publication_benchmark/current.json)
now admits 39/39 sourced potentials from 30 exact free molecules across seven
studies (MAE 0.247 V, RMSE 0.279 V). All outliers remain. This quantifies free-parent
accuracy and does not validate exact grafted structures or polymers.

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
| Oligomer graphs, occupancy and stereochemistry | [results/polymer_models/current.csv](polymer_models/current.csv); checked quantum outputs: [interaction status](polymer_models/interaction_status.json) and [spacing pointer](polymer_models/spacing_current.json) |

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

The [controlled CPU concurrency experiment](validation/cpu_concurrency_2026-10-22/README.md) checks 15 real descriptor executions with one to eight workers, including a two-host condition. All 12,555 molecular comparisons match serial definitions; launch overhead limits the speedup of this small library. The [shared route-stage launcher](validation/route_launcher_2026-10-22/README.md) also passes its 837-design execution. A [small CPU reproduction example](../examples/publication_smoke/README.md) is available, with [clean-checkout and fresh-environment reproduction checked](validation/fresh_checkout_2026-10-22/README.md).

The colleague-review revision has a [comment-by-comment ledger](../manuscript/revision_plan.md), a [rebuilt main paper](../manuscript/manuscript.pdf) and [separate SI](../manuscript/supporting_information.pdf). New checks include [grouped bias/reference diagnostics](validation/reference_bias_2026-10-22/README.md), [relative-basicity controls](validation/proton_transfer_2026-10-22/README.md) and [actual singlet/doublet quantum restart equivalence](validation/quantum_restart_2026-10-22/README.md). [Broader conformer pilots](validation/conformer_revision_2026-10-22/README.md) are active; DFT re-ranking and ensemble G remain pending. The [exact five-repeat phthalimide drawing](figures/polymer_models/spacing_construction.pdf) corresponds to the occupancy cycles. These checks do not make the project submission-ready.

A [standard-library cached-state example](../examples/cached_thermochemistry/README.md) now reconstructs one PMDI potential and two grafting sums after relocation, with no original NFS reads or scientific-package imports. [Recorded conda artifacts and pip versions](../environments/README.md) improve environment reproducibility, while full path portability, fresh GPU restoration and a public DOI remain open. All [five neutral-SVPD sensitivity controls](validation/neutral_basis_2026-10-22/README.md) are independently checked; use the independent review for completed quantities and leave pending fields unresolved.

The [all-axis decision publication](selection/README.md) now follows normal candidate finalization and figure updates. The fresh candidate reconstruction has 628 independent checks with no discrepancies. A [descriptor Pareto figure](figures/selection/descriptor_pareto.pdf) shows eight joint descriptor-front models; the [decision matrix](figures/selection/candidate_axis_matrix.pdf) separates supporting results from full-axis acceptance. No candidate is yet qualified on all required axes. The [targeted expansion plan](../docs/TARGETED_CANDIDATE_EXPANSION.md) defines 34 comparison slots, not newly enumerated molecules. [Twelve expanded CREST searches](validation/conformer_sampling_2026-10-23/README.md) have actual metadynamics starts; DFT and ensemble follow-up remain pending.

The [all-pose native baseline](polymer_models/campaigns/2026-10-23_baseline_two_pose_ensemble_v2/README.md) now has twelve checked G values and six independently reconstructed finite two-pose partitions. Adjacent/one-spacer nonadditivities are −5.17/−11.39 kJ/mol at the primary rotor cutoff; the [baseline figure](figures/polymer_models/conditional_ensemble_baseline.pdf) shows the prescribed cutoff sensitivity. [Broader search-to-DFT continuation](polymer_models/campaigns/2026-10-23_crest_refinement/README.md) remains active. Two older warning-heavy pilots were retired with all checkpoints preserved. Finite two-pose completion does not establish adequate conformer sampling.

The [exact PMDI candidate follow-up](validation/pmdi_occupancy_2026-10-23/README.md) links one enumerated design to reviewed precursor access, three checked polymer graphs, active source searches and a separate initial DFT baseline. The [construction drawing](figures/polymer_models/pmdi_occupancy_construction.pdf) shows the actual open graph. [Composition capacity bounds](figures/polymer_models/pmdi_capacity_bounds.pdf) include scaffold and residual-chloride mass; they are not measured accessible capacities or inferred loading distributions. These bounds enter the all-axis ledger as supporting scenario evidence, without a capacity pass or synthesis admission.

The [segment-deformation audit](reviews/2026-10-23/segment_deformation_audit_v2/README.md) reconstructs all three small totals and checks actual relaxation. The [individual-term figure](figures/polymer_models/segment_deformation_components.pdf) shows opposing fragment contributions in two pairs. These are signed differences from a sampled monomer reference, rather than positive-definite strain penalties.
