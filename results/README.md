# Results index

Use these stable paths. Completed quantities carry field-specific QC; blank values
are pending or excluded. Raw calculations and earlier source snapshots are retained.

The [current polymer execution update](validation/pcm_switch_curvature_guarded_2026-10-23/README.md) separates passed response controls from pending native G. Both exact 115-atom response tests and independent replays pass; no PMDI native G is yet accepted. The two-graft worker was stopped after a confirmed shared-device audit, with its 23 evaluated checkpoints preserved. A subsequent delayed foreign attachment was caught by the live guard; a fresh continuation is registered for an exclusively available suitable device. [Finite residency guards](validation/gpu_exclusivity_2026-10-23/README.md) now protect subsequent workers. These execution results do not admit a full-axis candidate.

| Result | Main entry | Current data |
| --- | --- | --- |
| Grafting thermochemistry: 75 reaction models | [Visual report and structures](grafting/index.html) | [Method/records](grafting/README.md) · [Status](grafting/status.json) · [Table](grafting/current.csv) |
| Polymer-fragment and interaction models | [3D interaction report](polymer_models/index.html) | [Method/records](polymer_models/README.md) · [Interaction status](polymer_models/interaction_status.json) · [Graph inventory](polymer_models/current.csv) |
| Graft spacing and conditional grafting framework | [Exact occupancy cases and report](polymer_models/campaigns/2026-10-21_spacing_cds_recovery/index.html) | [Reproduction guide](../docs/POLYMER_SPACING.md) · [Current pointer](polymer_models/spacing_current.json) · [Independent validation](polymer_models/workflow_validation.json) |
| Candidate redox properties and comparisons | [Candidates](candidates/README.md) | [Workbook](candidates/latest.xlsx) · [All candidates](candidates/figures/data/all_candidates.csv) |
| Maintained figures and molecular structures | [Figure guide](candidates/figures/README.md) | [Comparison browser](candidates/figures/index.html) |
| Measured/calculated redox validation | [Validation guide](validation/README.md) | [Figure](figures/validation/matched_parent_validation.pdf) · [Actual values](figures/validation/matched_parent_validation.csv) |
| Publication preparation and benchmark completion | [Readiness guide](publication_readiness/README.md) | [Requirements](publication_readiness/current.json) · [Frozen benchmark](validation/publication_benchmark/README.md) |
| Candidate search, selection and source datasets | [Design/evaluation report](discovery/expansion/index.html) | [Discovery guide](discovery/README.md) · [Matched effects](discovery/expansion/matched_design_effects.pdf) · [Search inventory](discovery/systematic_families_2026-10-18/summary.json) |
| Full property-table arithmetic review | [Reviewed tables](validated/README.md) | [Dated review pointer](validated/current.json); later candidates use the active candidate feed |
| Working in silico manuscript | [Draft PDF](../manuscript/manuscript.pdf) | [Source and build](../manuscript/README.md) · [Early review questions](../manuscript/review_notes.md) · [Frozen evidence](../manuscript/evidence/2026-10-23_thermal_cutoff_revision/manifest.json) |
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


The October 23 manuscript revision adds conjunctive admission across all required
axes, an exact PMDI precursor/route/capacity assessment and all twelve inherited
phthalimide native free energies in a finite two-pose partition. The checked
paper has 15 pages, including its TOC page, and the SI 12. The [compact PMDI construction](figures/polymer_models/pmdi_occupancy_construction_compact.pdf) preserves the
prepared graph while making chemical labels legible at manuscript size. No full
material candidate is admitted; broader sampling and exact PMDI occupancy
calculations remain separate work.


The [portable scientific-result package](releases/primary_review/README.md) now
contains 243 thermochemical states and 32 electronic-only records. It replays
primary and auxiliary SI quantities after relocation; 242 Hessian partitions
and both restart geometry comparisons pass CPU reconstruction. All 2,498
archive files are verified. Public deposition, diagnostic artifact closure and
full GPU restoration remain open.

The [DIPEA reference-origin audit](validation/dipea_reference_origin_2026-10-23_checked/README.md) excludes the reported MeCN pKaH 18.6 from independent measured validation because its specific source note identifies a PCM estimate while another table gives a generic experimental label. The current manuscript records this conflict. Primary reaction energies and candidate decisions are unchanged.

The [32 GB solvent-workspace qualification](validation/smd_workspace_32GB_2026-10-23/README.md) records passed small spin controls, completed 115-atom single-point SCF checkpoints and the subsequent preserved gradient allocation failure. The separate bounded-gradient control also passes. Large-geometry optimization and native thermochemistry on this hardware remain unqualified.

The [large 115-atom gradient recovery](validation/smd_gradient_workspace_32GB_2026-10-23/README.md) now passes four sampled total-energy derivative checks on a 32 GB V100, with 20,304 MiB peak GPU memory. Its same-model [full geometry/energy continuation](polymer_models/campaigns/2026-10-23_crest_refinement/baseline__design_pmdi_ethyl__occ_2_workspace_32GB/README.md) has an actual reserved optimizer start; native G and ensemble convergence remain pending. The candidate-only observer resumed after the documented handoff; active claims suppress duplicate workers. The independent wider-search observer and source searches continue.

A separate [chemical TOC graphic](figures/publication/README.md) is prepared from exact precursor/product identities, with EPS and 600 dpi color TIFF exports. It is visually checked in color and grayscale and inserted on the checked manuscript's final page; author review remains pending.

The [current manuscript revision](../manuscript/README.md) adds explicit reviewer-package scope and the integrated TOC. Its frozen counts are independently cross-checked against the archived inventory. The [site-3 post-pilot continuation](polymer_models/campaigns/2026-10-23_crest_refinement/baseline__design_pmdi_ethyl__occ_3_postpilot/README.md) resolves a force-directory collision without moving the original pilot. The [64 GB profile campaign](validation/smd_gradient_workspace_64GB_2026-10-23/README.md) has passed both small spin controls and started the separate 139-atom resource requalification. Native polymer G remains pending.

The original [two-graft 139-atom full continuation](polymer_models/campaigns/2026-10-23_crest_refinement/baseline__design_pmdi_ethyl__occ_2_3_workspace_64GB/README.md) passed a 39,206 MiB point qualification but later stopped at an XC-gradient workspace allocation. Its raw trajectory is preserved. The [workspace execution workflow](../docs/POLYMER_WORKSPACE_EXECUTION.md) supports fresh full-calculation storage and checked-energy handoff to native thermochemistry. Site 3 now has checked electronic energy and a separate active native-thermal pilot; site 2 stopped after failing the final total-force energy-difference check. G remains pending until the normal thermal checks pass.

The checked v5 archive is now uploaded as an [unpublished draft release](releases/primary_review/draft_release.json). Remote downloads of all three assets match local checksums. It precedes current polymer calculations and does not constitute public deposition or a DOI.

The [rotor-cutoff sensitivity audit](validation/thermal_cutoff_sensitivity_2026-10-23/README.md) independently reconstructs all 39 benchmark endpoints and 87 completed grafting reaction models at three fixed-geometry settings. Its [two-panel figure](figures/validation/thermal_cutoff_sensitivity.pdf) is included in the current SI. Maximum spans are 25.43 mV and 4.46 kJ/mol, with no nominal reaction sign changes; no missing chemical criterion is inferred to pass.

The [XC temporary-lifetime controls and checkpoint check](validation/xc_workspace_lifetime_2026-10-23/README.md) pass both small spin controls and sampled total-energy derivatives at the preserved 139-atom checkpoint. A [fresh full continuation](polymer_models/campaigns/2026-10-23_crest_refinement/baseline__design_pmdi_ethyl__occ_2_3_xc_checkpoint_resume/README.md) is registered for the H100 node; its actual worker start remains separate from preparation. The [site-3 native-thermal pilot](polymer_models/campaigns/2026-10-23_crest_refinement/baseline__design_pmdi_ethyl__occ_3_native_xc_surface_budget/README.md) is running on a V100 using the measured-cavity host forecast and unchanged scientific checks. Neither branch supplies a qualified G yet.

The [site-2 PMDI surface audit](validation/pmdi_site2_CDS_surface_2026-10-23/README.md)
quantitatively attributes its failed coarse-step force comparison to the legacy
SMD CDS term. Independent replay checks 29 energy/coordinate samples; the
[diagnostic figure](figures/validation/pmdi_site2_CDS_surface.pdf) separates this
term from the remaining discrepancy. A [separate full quantum fine-step test](validation/pmdi_site2_fine_force_2026-10-23/README.md)
has an actual reserved V100 start and durable evaluations. The original state
remains unaccepted; no native G, conformer convergence or all-axis pass follows
from the component audit.

The site-2 full quantum fine-step diagnostic has now passed all 15 evaluations and independent replay. A [separate electronic qualification](polymer_models/campaigns/2026-10-23_crest_refinement/baseline__design_pmdi_ethyl__occ_2_fine_force_qualified/README.md) records the named four-step validation policy and unchanged tolerance while preserving the original failed workflow. Its saved E can now enter the ordinary native-thermal handoff; no native G, ensemble convergence or synthesis admission is inferred.

The site-3 V100 native pilot has now failed in the analytic PCM Hessian response, with its raw records preserved and no accepted G. A [fresh H100 native qualification](polymer_models/campaigns/2026-10-23_crest_refinement/baseline__design_pmdi_ethyl__occ_3_native_xc_H100/README.md) is registered through the existing observer. The separate site-2 native V100 pilot remains active; no large-native hardware qualification is inferred from small spin or fixed-gradient controls.

The [host switch-curvature qualification](validation/pcm_switch_curvature_guarded_2026-10-23/README.md) addresses a later 43.37 GiB allocation that response-only tests do not cover. Both guarded full spin/Hessian controls pass. Fresh native pilots retain all physical gates and an explicit measured-peak/device-headroom policy; their G values are still pending. The [context revision](../manuscript/README.md) adds a checked aqueous NDI association reference without changing the MeCN evidence epoch.
