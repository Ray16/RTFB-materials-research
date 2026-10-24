# Results index

Use these stable paths. Completed quantities carry field-specific QC; blank values
are pending or excluded. Raw calculations and earlier source snapshots are retained.

The independent exact-PMDI reviews verify both [site-2](reviews/2026-10-24/pmdi_site2_native/README.md) and [site-3](reviews/2026-10-24/pmdi_site3_native/README.md) native free-energy components. The 139-atom two-graft model has checked stationary coordinates and canonical optimized-state gas/SMD energies. Its native free energy, matched occupancy cycle and broader conformer convergence remain pending; no full-axis candidate is admitted.

The [139-atom V100 reference point](polymer_models/campaigns/2026-10-23_crest_refinement/baseline__design_pmdi_ethyl__occ_2_3_host_workspace_32GB/README.md), recovered final geometry checks and canonical optimized-state gas/SMD energies now pass. The [native-response probe](validation/pmdi_pair_PCM_response_v2_2026-10-24/failure_review.json) preserved its converged baseline after a large contraction exceeded GPU memory. The [tiled FP64 contraction implementation](validation/streamed_weighted_2026-10-24/README.md) passes both full spin/Hessian controls; its exact-model response and derivative check also passes independent reconstruction. The [full paired native calculation](polymer_models/campaigns/2026-10-24_solute_workspace_native/baseline__design_pmdi_ethyl__occ_2_3/README.md) is running after additional [pure-solute workspace controls](validation/solute_workspace_2026-10-24/README.md); G and the matched cycle remain pending.

The [current UMA entry point](validation/uma_current_2026-10-24/README.md) now prefers Nitro and has a tested V100 fallback. The [site-2 PMDI native component](reviews/2026-10-24/pmdi_site2_native/README.md) independently reconstructs, bringing the singly occupied PMDI components to two; the matched two-graft cycle remains pending.

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
| Working in silico manuscript | [Draft PDF](../manuscript/manuscript.pdf) | [Source and build](../manuscript/README.md) · [Early review questions](../manuscript/review_notes.md) · [Frozen evidence](../manuscript/evidence/2026-10-24_spatial_revision/manifest.json) |
| Research presentation | [Slides](../slides/README.md) | [Checked PDF](../slides/redox_workflow.pdf) |
| Detailed validation/provenance index | [Latest review](LATEST_REVIEW.json) | [Machine-readable result index](INDEX.json) |
| V100 inference compatibility and loader behavior | [V100 guide](validation/uma_v100_2026-10-24/README.md) | [Current status](validation/uma_v100_2026-10-24/summary.json) · [Unaltered comparison](validation/uma_v100_2026-10-24/reference_comparison.json) · [Composition diagnostic](validation/uma_v100_2026-10-24/composition_comparison.json) |
| Streamed V100 solvent-Hessian execution | [Method and running polymer recovery](validation/streamed_PCM_2026-10-24/README.md) | [Full spin replays](validation/streamed_PCM_2026-10-24/replay_qualification.json) · [Large-case start](polymer_models/campaigns/2026-10-24_streamed_native/baseline__design_pmdi_ethyl__occ_2/startup_verification.json) |
| Accepted spacing-control conformers | [Actual molecular geometries](figures/polymer_models/spacing_conformers.pdf) | [Coordinates, contact metrics and source hashes](figures/polymer_models/spacing_conformers.json) |

The [PMDI attachment-class review](reviews/2026-10-24/pmdi_alkylation_context/README.md) adds primary molecular N-alkylation context. It preserves the exact-route hold and all 75 route outcomes.

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

A [standard-library cached-state example](../examples/cached_thermochemistry/README.md) now reconstructs one PMDI potential and two grafting sums after relocation, with no original NFS reads or scientific-package imports. [Recorded conda artifacts and pip versions](../environments/README.md) improve environment reproducibility, while full path portability and a public DOI remain open; the later scoped fresh-GPU controls are indexed below. All [five neutral-SVPD sensitivity controls](validation/neutral_basis_2026-10-22/README.md) are independently checked; use the independent review for completed quantities and leave pending fields unresolved.

The [all-axis decision publication](selection/README.md) now follows normal candidate finalization and figure updates. The fresh candidate reconstruction has 628 independent checks with no discrepancies. A [descriptor Pareto figure](figures/selection/descriptor_pareto.pdf) shows eight joint descriptor-front models; the [decision matrix](figures/selection/candidate_axis_matrix.pdf) separates supporting results from full-axis acceptance. No candidate is yet qualified on all required axes. The [targeted expansion plan](../docs/TARGETED_CANDIDATE_EXPANSION.md) defines 34 comparison slots, not newly enumerated molecules. [Twelve expanded CREST searches](validation/conformer_sampling_2026-10-23/README.md) have actual metadynamics starts; DFT and ensemble follow-up remain pending.

The [all-pose native baseline](polymer_models/campaigns/2026-10-23_baseline_two_pose_ensemble_v2/README.md) now has twelve checked G values and six independently reconstructed finite two-pose partitions. Adjacent/one-spacer nonadditivities are −5.17/−11.39 kJ/mol at the primary rotor cutoff; the [baseline figure](figures/polymer_models/conditional_ensemble_baseline.pdf) shows the prescribed cutoff sensitivity. [Broader search-to-DFT continuation](polymer_models/campaigns/2026-10-23_crest_refinement/README.md) remains active. Two older warning-heavy pilots were retired with all checkpoints preserved. Finite two-pose completion does not establish adequate conformer sampling.

The [actual finite-control conformers](figures/polymer_models/spacing_conformers.pdf) show close and separated cores for both sequence spacings. Within those two-pose sets, close-contact conformers have lower G by 11.9 and 7.1 kJ/mol for adjacent and one-spacer graphs. These differences include the entire chain and solvent; they do not isolate motif binding or prove converged populations. The plotted backbone follows the exact checked XYZ connectivity and uses proper rotations without reflection or coordinate rescaling.

The [exact PMDI candidate follow-up](validation/pmdi_occupancy_2026-10-23/README.md) links one enumerated design to reviewed precursor access, three checked polymer graphs, active source searches and a separate initial DFT baseline. The [construction drawing](figures/polymer_models/pmdi_occupancy_construction.pdf) shows the actual open graph. [Composition capacity bounds](figures/polymer_models/pmdi_capacity_bounds.pdf) include scaffold and residual-chloride mass; they are not measured accessible capacities or inferred loading distributions. These bounds enter the all-axis ledger as supporting scenario evidence, without a capacity pass or synthesis admission.

The [segment-deformation audit](reviews/2026-10-23/segment_deformation_audit_v2/README.md) reconstructs all three small totals and checks actual relaxation. The [individual-term figure](figures/polymer_models/segment_deformation_components.pdf) shows opposing fragment contributions in two pairs. These are signed differences from a sampled monomer reference, rather than positive-definite strain penalties.


The current manuscript includes conjunctive admission across all required
axes, an exact PMDI precursor/route/capacity assessment and all twelve inherited
phthalimide native free energies in a finite two-pose partition. The checked
paper has 17 pages, including its TOC page, and the SI 14. The [compact PMDI construction](figures/polymer_models/pmdi_occupancy_construction_compact.pdf) preserves the
prepared graph while making chemical labels legible at manuscript size. No full
material candidate is admitted; broader sampling and exact PMDI occupancy
calculations remain separate work.


The [portable scientific-result package](releases/primary_review/README.md) now
contains 243 thermochemical states and 32 electronic-only records. It replays
primary and auxiliary SI quantities after relocation; 242 Hessian partitions
and both restart geometry comparisons pass CPU reconstruction. All 2,498
archive files are verified. Public deposition, diagnostic artifact closure and
broader deployment portability remain open; scoped fresh-GPU controls are indexed below.

The [DIPEA reference-origin audit](validation/dipea_reference_origin_2026-10-23_checked/README.md) excludes the reported MeCN pKaH 18.6 from independent measured validation because its specific source note identifies a PCM estimate while another table gives a generic experimental label. The current manuscript records this conflict. Primary reaction energies and candidate decisions are unchanged.

The [solvent-workspace controls](validation/pcm_switch_curvature_guarded_2026-10-23/README.md), [site-2 fine-step force qualification](polymer_models/campaigns/2026-10-23_crest_refinement/baseline__design_pmdi_ethyl__occ_2_fine_force_qualified/README.md) and [execution workflow](../docs/POLYMER_WORKSPACE_EXECUTION.md) retain the numerical evidence and prior failed attempts. The current [site-3 free-energy review](reviews/2026-10-24/pmdi_site3_native/README.md) and [split geometry/energy branch](validation/pmdi_split_execution_2026-10-24/README.md) own the latest scientific and execution status. Historical control success alone does not supply a missing large-state energy or G.

A separate [chemical TOC graphic](figures/publication/README.md) is prepared from exact precursor/product identities, with EPS and 600 dpi color TIFF exports. It is visually checked in color and grayscale and inserted on the checked manuscript's final page; author review remains pending.

The checked v5 archive is now uploaded as an [unpublished draft release](releases/primary_review/draft_release.json). Remote downloads of all three assets match local checksums. It precedes current polymer calculations and does not constitute public deposition or a DOI.

The [rotor-cutoff sensitivity audit](validation/thermal_cutoff_sensitivity_2026-10-23/README.md) independently reconstructs all 39 benchmark endpoints and 87 completed grafting reaction models at three fixed-geometry settings. Its [two-panel figure](figures/validation/thermal_cutoff_sensitivity.pdf) is included in the current SI. Maximum spans are 25.43 mV and 4.46 kJ/mol, with no nominal reaction sign changes; no missing chemical criterion is inferred to pass.

The [fresh quantum deployment](validation/gpu_environment_restore_2026-10-23/README.md) now reproduces both complete spin controls with independent Hessian/partition replay. A separate compatibility profile fixes fourteen recorded CUDA dependency conflicts while retaining method-package versions and the xTB hash. The [manuscript](../manuscript/manuscript.pdf) integrates this checked result and frames the chemical findings without execution-history narration. Exact-candidate occupancy, broader conformer convergence and public deposition remain open.

The [UMA Nitro validation](validation/uma_nitro_2026-10-23/README.md) has a pinned separate environment and a ten-case exact-state panel. Legacy inference reproduces six fingerprinted stored energies. [New CUDA12 FP32 inference runs on V100](validation/uma_v100_2026-10-24/README.md); a controlled composition-averaging ablation explains its differences from legacy inference across all forty frames. The diagnostic override is not a production configuration. Actual Nitro energy/force and matched-device speed validation of the redox panel remain pending. The current public CLI prefers Nitro for supported molecular warm starts, uses the tested legacy V100 fallback and keeps the model profiles separate; see the [current execution policy](validation/uma_current_2026-10-24/README.md). [CPU electronic fallback controls](validation/CPU_quantum_2026-10-23/README.md) pass two spin comparisons at the unchanged quantum level; the large exact-polymer resource point remains separate.

The [first final-CREST geometry refinement on V100](polymer_models/campaigns/2026-10-24_CREST_geometry/alternating__occ_none__conf001/README.md) is running as a separate stage after an independently reviewed exact-input resource point. The low-level search contains 4,329 topology-checked conformers; twenty starts are staged. Neither DFT re-ranking nor a broader ensemble free energy is complete.

Both empty-backbone searches now have twenty staged starts per configuration. The [second configuration resource point](polymer_models/campaigns/2026-10-24_CREST_V100/all_same__occ_none_point80/README.md) is running on a separate V100. Complete DFT re-ranking and broader finite ensembles remain pending.
