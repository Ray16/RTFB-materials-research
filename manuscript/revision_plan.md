# Colleague review: revision ledger

Review received 22 October 2026, referring to the immutable 21 October draft.
The review's scientific concerns remain applicable even where counts have since
changed. A task is closed by checked evidence and a corresponding manuscript
revision, not by a proposed calculation or a running job.

The separate [TEMPO exact-support assessment](../results/grafting/tempo_exact_support_2026-10-26/README.md)
adds a source-supported attachment and independently reconstructed local-model
oxidation potential. Its route passes, but nine required material axes remain
pending. The checked TEMPO and local-motif raw-state addenda are included in
the locally replayed consolidated reviewer archive; public deposition remains open.

## Major comments and submission requirements

| Comment | Current evidence and required revision | Status |
| --- | --- | --- |
| M1: no end-to-end decision | Exact enumerated PMDI and NDI matches now have freshly checked six-state redox/Fc reconstruction, grafting hashes and two mapped route assessments (both hold). Exact PMDI precursor access is source-checked and the site-2 and site-3 polymer native G now pass independent reconstruction. Balanced first-site DBU reactions and polymer-minus-local penalties now independently reconstruct and appear in the paper/SI. The complete four-state PMDI native cycle now independently reconstructs at −1.86 kJ/mol, with the actual separated-core geometry and both conditional attachment orders reported. Broader conformer sampling, selectivity and core retention remain unresolved. Complete the missing assessments and present which evidence changes priority or leaves synthesis on hold. Method controls must not manufacture synthesis eligibility. | Finite exact PMDI worked example now reported; ensemble/material selection remains open |
| M2: benchmark bias/reference | Keep all 39 unfitted endpoints primary. Report the existing molecule-cluster intervals, qualify the seven unequal study clusters, and add endpoint-, molecule- and study-held-out offset analyses, including redox-step-specific models. Compare Fc-independent spacings and multiple organic-reference diagnostics. Verify the Fc literature's methods and conventions. Discuss the orbital/Jahn–Teller limitation without assigning an unvalidated correction. | Analysis and six fixed-geometry functional controls checked and included; reference attribution remains conditional |
| M3: constant base shift | The per-base constant and DBU/HCl discrepancy are reported. A 2021 spectrophotometric study reports DIPEA log Kᴴ = 18.64(1) in MeCN at 25 °C; a separate 2024 NMR study estimates 18.7 for DIPEAH hexafluorophosphate in MeCN-d3. The conditional cross-source DIPEA/DBU transfer comparison is +31.80 computed versus +32.36 kJ/mol implied by the 18.64 result; the DIPEA/HCl capture terms are −70.98 versus −47.60 kJ/mol. Indicator calibration, counterion and activity conventions are not shown to match, and HCl cancels from the complete DFT grafting cycle. Hybrid replacements are sensitivities, not validated full-reaction corrections. The older 18.6 PCM source stays separately identified. | Two experimental source estimates checked; full-route/speciation qualification open |
| M4: conformer/tacticity sensitivity | All twelve inherited phthalimide pose G values now pass independent reconstruction; the finite two-pose ensemble gives -5.17/-11.39 kJ/mol, but this does not establish sampling convergence. Document seed 20261020, 4 embeds, 2 xTB seeds and 2 DFT poses per graph. Add broader conformer search, distinct backbone stereochemical scenarios, DFT re-ranking and a defined ensemble treatment. Predeclare convergence criteria. Until demonstrated, retain the result as a finite method-control example, not a spacing design rule. | Open; highest scientific priority |
| M5: geometry/Hessian basis imbalance | Five matched neutral-SVPD controls independently reconstruct, giving -11.7 to -35.1 mV shifts with the existing anion/Fc fixed. The main paper and SI now report geometry, nuclear-partition and resulting potential changes separately. Keep primary source-frozen records unchanged. Alternative-basis controls must pass their physical checks and cannot inherit the original certificate by relabeling. | Five-control sensitivity checked and reported; full basis convergence not claimed |
| M6: Methods/QC precision | State exactly two tested Hessian columns and both displacement sizes; explain |nu| for tolerated -5..0 cm^-1 modes, unscaled frequencies, spin-contamination threshold 0.1, actual force values and xTB symmetry detection. Remove the unsupported adjective “conservative” for the 4.5e-4 force threshold. Cite qRRHO and program/method papers. | Exact gate description included; further physics controls remain separate |
| M7: portable/public reproduction | The local reviewer bundle includes 243 thermochemical states and 32 electronic-only records. Relocated offline arithmetic, 242 Hessian reconstructions and seven negative checks pass; the current snapshot freezes and reconciles these inventories. Public code visibility is verified and environment artifacts are pinned. A separate fresh compatible GPU environment now reproduces two complete spin controls and independent Hessian partitions. Public raw-state deposition/DOI, final metadata, remaining diagnostic closure and broader path portability remain open. | Locally prepared and checked; public/fresh-deployment requirements open |
| M8: literature context | Verify each DOI and source. Add computational electrochemistry/reference, high-throughput design, redox-polymer, workflow/provenance and missing method literature where it supports a specific argument. Do not turn remembered references into unverified citations. | Bibliography now has 36 entries, including the exact precursor patent and both primary DIPEA experimental sources; further relevant context in progress |
| M9: cohesion and hedging | Lead with the chemical decision enabled by the identity-linked framework. Move synthetic scheduling, concurrency, restart and numerical-repair detail to a separately built SI. Consolidate limitations, define unfamiliar terms and reduce implementation vocabulary in the main article. | Main and separate SI built; 19 (including the TOC page) and 20 pages, respectively |
| Remaining endpoint | The fixed panel now has 39/39 endpoints; the last triplet alternative is 0.707 eV above the selected singlet electronically. Independent reconstruction: 111 checks, zero discrepancies. Alternative-spin thermochemistry remains explicitly unresolved. | Calculation/reconstruction resolved; latest draft updated |

An [external NDI measurement](../results/validation/ndi_external_challenge_2026-10-26/README.md)
has been source-checked in MeCN, with an exact proposed molecular graph and
explicit quasi-reversible peak limitations. Its states have not been computed;
it is not in the fixed benchmark and does not yet establish NDI accuracy.
Ten graph-checked [force-field starting geometries](../results/validation/ndi_external_control_seeds_2026-10-26/README.md)
are prepared for the exact source molecule; the state-specific calculations
remain pending.
An isolated [UMA control adapter](../results/validation/ndi_external_uma_2026-10-26_v2/README.md)
is now prepared and tested without a GPU. A neutral-state resource pilot is
gated on independent fixed-device Nitro parity; all state-specific results and
the broader dispatch remain pending.

The broader-conformer continuation now has declared refinement coverage for all
240 separated starts across twelve graph/configuration sets. The last two
occupied geometry resource classes passed their independent reviews. Separate
occupied hybrid-energy pilots and an empty-chain native pilot are still physical
validation steps, not finished ensembles. A source-bound continuation requires
the native pilot's independent review and each exact electronic parent before
launching additional native states; its dependency/resource checks and the native
bookkeeping pass 38 focused tests. This leaves M4 open and does not change the
checked manuscript's finite-sample claims. See
[the current follow-up checkpoint](../results/reviews/2026-10-26/conformer_followup_milestone/README.md).

The separate smooth-surface 91-atom native pilot now passes independent full
Hessian/mode/nuclear reconstruction, enabling the forty-state continuation.
The occupied 105-atom hybrid point also passes independent checkpoint/resource
review and feeds a separate native pilot. Neither result closes ensemble
convergence. A new frozen methyl-capped phthalimide control supplies the first
direct local interaction triplet (−24.76 kJ/mol in MeCN; −37.16 in gas).
All twelve pose/cap interaction triplets now reconstruct from 36 independently
reviewed electronic states; the second close pose is −26.22 kJ/mol in MeCN at
the central cap length, and the separated poses are near zero. Basis-superposition
checks and broader conformer sampling remain pending.
This metric includes caps and continuum response and remains distinct from
whole-chain ΔΔG, pure π binding and electron-transfer coupling. The current
manuscript and SI now report the complete frozen local panel with these limits. See
[the local interaction record](../results/reviews/2026-10-26/local_motif_interactions/README.md).
Four checked empty-chain native records in the alternate smooth-surface model
now span 7.19 kJ/mol at the 50 cm^-1 rotor cutoff. This is a geometry sensitivity
within one configuration, not a completed ensemble, grafting-cycle uncertainty
or a result to mix with the primary-model cycle. See the [range diagnostic](../results/validation/empty_chain_native_range_2026-10-26/README.md).
The four geometries pass graph/stereochemistry checks and are separated by
2.01--4.27 Å heavy-atom RMSD; a finite-set partition lowers the minimum by
1.81 kJ/mol. This does not satisfy the 10-versus-20 sampling criterion. See the
[finite-set analysis](../results/validation/empty_chain_finite_set_2026-10-26/README.md).
The occupied site-2 native pilot was independently reviewed, and its
[follow-up coordinators](../results/reviews/2026-10-26/occupied_followup_resumption/README.md)
have been resumed with dependency and reservation checks. Twelve prior intents
were subsequently verified as unstarted and reopened. These are process and input records, not
completed ensemble data.
The electronic adapter's remote working-directory failure has since been
reproduced and corrected; see the [placement recovery](../results/reviews/2026-10-26/electronic_workdir_recovery/README.md).
The [complete local-panel review](../results/reviews/2026-10-26/local_motif_complete/README.md)
preserves the earlier partial snapshot and the full checked replay.

The [paired-cycle raw addendum](../results/releases/neighbor_cycle/README.md)
closes the latest local archive-coverage gap: fourteen cached states and balanced
comparisons reproduce the complete exact-PMDI cycle at all three rotor cutoffs
after extraction with original checkout reads blocked. Seven nominal/adverse-input
checks initially passed. The extended package now also reconstructs thirteen
molecular Hessian mode sets and nuclear sums at three cutoffs, plus monatomic
chloride translation, in the locked fresh NumPy environment. Fourteen nominal
and adverse-input checks pass. This extends the earlier thirteen-state first-site
artifact without altering it. Public deposition, general cluster-path portability
and conformer convergence remain open. The revised manuscript and SI report the
new reproducibility evidence without changing the chemical numbers or scope.

The end-to-end examples and adequate conformer evidence are required before
claiming a complete selection framework. The open PMDI trimer slots and an NDI
external anchor need not become headline claims. Existing calculations remain
preserved and their scientific limits remain visible.

## Minor comments and exact paper/code discrepancies

| Comment | Required disposition |
| --- | --- |
| m1: stale draft | Latest draft uses a new immutable snapshot; 39/39 benchmark and 6/6 selected spacing G. Preserve old snapshots rather than editing their evidence. |
| m2: cluster intervals | Put intervals in the paper and state that seven unequal study clusters provide weak interval calibration. Report molecule and study resampling separately. |
| m3: wrong model in construction figure | Draw the exact phthalimide five-repeat occupancy model, including styrene padding, remaining benzyl chlorides, mapped sites and stereochemistry. |
| m4: two 75-member cohorts | Join the candidate and reaction ledgers explicitly; distinguish 75 candidate reactions from 72 unique precursors and explain exact identity overlap. |
| m5: small segment deformation | Report the numbers in text/SI and audit what relaxation and reference geometries actually occurred; do not infer a relaxed ensemble from a small deformation value. |
| m6: conversion constant | Use exact SI/CODATA conversion in new analysis; preserve and label the frozen 96.485-source arithmetic rather than silently changing old values. |
| m7: build provenance | Separate compiler-produced checks from a dated, hash-bound human visual-review record. Do not imply manually added review fields were emitted by the build. |
| m8: figures/graphical abstract | Chemical TOC artwork uses exact precursor/product graphs and a schematic neighbor panel; inserted on the final manuscript page. EPS/TIFF dimensions, independent rendering and final PDF layout are checked. Author approval remains pending. |
| m9: jargon | Replace implementation labels where possible and define the remaining evidence/decision terms in one compact table. |
| m10: reproducibility details | Record seed, sampling counts, backbone configuration and xTB symmetry detection explicitly. |
| Hessian replay rounding | Correct 1.34e-11 to the value supported by the stored evidence (1.333e-11; 1.33e-11 at three significant figures). |
| QC/imaginary-mode treatment | Match the text to actual sampled columns, thresholds and tolerated-negative-frequency handling. |
| Base-induced sign changes | Show the per-base shift across the whole qualified panel, including AQ and NDI, without a ranking-change claim. |
| xTB build/path locking | Pin and document the actual quantum environment; configure deployment paths without silently bypassing scientific integrity or admission. |

## Work already running when this review arrived

The stacked-trimer optimizer returned convergence and an active SP, then a later library call crashed. Separate fixed/fine final-force checks failed, so the case remains excluded. Its obsolete finite-panel monitors have stopped; the two broader CREST searches continue independently. Its small-step force diagnostic is complete. They remain distinct from the completed six-state spacing
panel. The two-spin quantum interruption experiment has completed and passed an
independent reconstruction of its four final thermal partitions. The four relative-basicity components and both balanced transfer sums have passed independent reconstruction. The two original low-thread empty-backbone pilots have been retired with their raw checkpoints preserved. Of twelve broader searches, four completed and eight reached their declared timeout; all source processes ended. A new graph/stereochemistry-aware selection supplies twenty separated input geometries for each of twelve graph/configuration sets. A separately controlled smooth-area model passes four complete small quantum controls and a 91-atom geometry pilot; forty new empty-backbone geometries are prepared with actual starts verified. Distinct models must not be mixed in a reaction cycle. Broader DFT re-ranking and converged ensemble G remain pending.
These execution checks belong in the SI and do not substitute for M1/M4.

## All-axis requirement and focused expansion, 23 October

The user's clarified rule is conjunctive acceptance on every required axis before
Pareto comparison. The new source-current delivery stage implements this with
identity, model-scope, context and interval checks. Supporting proxies do not
supply missing properties. The current profile still has four undeclared
application cutoffs and incomplete physical evidence, so 878 candidate designs
remain pending and one comparison control is separate. The eight-model molecular
descriptor front is explicitly not a full material shortlist. The 34-slot focused
expansion plan is predeclared comparison coverage, not 34 additional molecules.
These improvements do not close the end-to-end or polymer-convergence comments.

A separate exact AQ-2-carboxylate capped-site proposal has now been constructed
and route-screened against a soluble-polymer literature precedent. It remains
on hold with no new redox calculations, and does not complete the matched
34-slot plan or establish an all-axis decision. The source-linked graph and
gate are recorded in [the proposal](../results/discovery/aq_carboxylate_2026-10-26/README.md).

## Conformer thermochemistry follow-up, 23 October

A finite partition implementation and an all-pose native thermal adapter now
reconstruct the six inherited spacing G records against their active electronic
and projected nuclear/spin parents. All twelve original pose free energies are now complete and independently
reconstructed. The finite two-pose ensemble gives adjacent/one-spacer
nonadditivity of -5.17/-11.39 kJ/mol at the primary rotor cutoff. Prescribed
25/50/100 cm^-1 cutoffs give -5.13/-5.17/-4.64 and -12.34/-11.39/-11.69 kJ/mol,
respectively. This quantifies the inherited finite baseline while wider CREST
searches run; it does not bound the error from undiscovered conformers.
All six original pairs have heavy-atom RMSDs above 3 Å. Duplicate grouping uses
complete-link geometry/energy checks with cutoff sensitivity; search hit counts
never become statistical weights. Neither two-pose sampling nor a running job
closes M4. The main manuscript and its immutable evidence remain unchanged until
new physical results are independently checked.

The read-only segment-deformation audit resolves the numerical part of m5. The
small +0.38–0.56 kJ/mol totals are signed reference differences, and two include
cancellation between negative and positive fragment contributions. Pair relaxation
is 0.25–0.70 Å over 57–93 logged optimizer steps. The next SI revision must report
the individual terms and avoid calling the total a positive-definite strain
penalty. This remains a sampled, two-level electronic decomposition, with no
global monomer minimum or Hessian-minimum claim.


## Exact PMDI follow-through and material-capacity bounds, 23 October

The mono-N-ethyl/NH PMDI precursor now has an exact primary-source identity
check (EP0010440B1, Example 19, printed page 10). This closes precursor access
for this design only. Attachment selectivity and redox-core retention remain
unresolved, so the prospective route remains on hold and full synthesis
admission is false. Three exact five-repeat PMDI occupancy graphs have passed
connectivity, stereochemistry and starting-geometry review, and entered CREST
searches. Separate one-minimum DFT resource pilots are running; no PMDI
occupancy free energy is yet qualified. The largest initial pilot exhausted
GPU memory; its failed record is preserved. A temporary-buffer allocation change
passed original-versus-modified energy, geometry, gradient and Hessian controls
for both spin classes. The same-model larger-system recovery remains a resource
pilot, not a completed scientific result.

A composition model now independently matches the masses of all eight exact
occupancy graphs. At three reactive handles per five repeats, one/two/three
grafts give two-electron capacity ceilings of 61.2/98.9/124.5 mAh/g under the
declared neutral linear-fragment convention. These are conditional upper bounds,
not measured capacity or evidence of a complete redox path in the polymer.
Actual loading, accessible utilization and the application target remain
unresolved. The all-axis ledger retains these as supporting observations,
without converting the capacity axis or any candidate into a pass.

The final twelve-state audit, exact construction and capacity figures, updated
72-page research deck and software checks are maintained with the results.
The primary manuscript and immutable evidence have not been silently replaced.
M1 and M4 remain open until the corresponding physical evidence and manuscript
revision are complete.


## Manuscript integration, 23 October

A new immutable epoch (`2026-10-23_all_axes_checked`) incorporates the twelve
finite-partition records and independent review, exact PMDI precursor source,
composition bounds and all-axis ledger. The 14-page main paper now explains the
limited three-core coverage, conjunctive admission and the exact PMDI pending
decision; the 10-page SI reports the finite partition and signed fragment
deformation audit. M1 remains incomplete without accepted exact-candidate
occupancy results and the remaining chemical evidence. M4 remains incomplete
without wider conformer/configuration convergence. The final graphical abstract,
public raw-state deposition and author declarations remain open. No submission
readiness is inferred from a successful build.


## Portable primary raw-state package, 23 October

A source-bound 234-state package now reproduces 39 benchmark endpoints, 67
reference reactions, 20 base-assisted values, two relative-basicity controls
and the twelve-state finite polymer partition after relocation. The standard
library replay runs from `/` with original checkout/calculation opens blocked;
four corrupt/misbound/path-escape cases are rejected. An optional CPU check
reconstructs 233 Hessian mode sets and nuclear partitions. All 2,234 archived
files match their source hashes. This substantially advances M7 but does not
close the auxiliary SI raw-state coverage, fresh GPU deployment, data license,
author metadata or public DOI deposition. The manuscript's immutable evidence
and checked PDFs are unchanged during this release-preparation pass.


## Auxiliary scientific records and memory recovery, 23 October

The review package now includes five alternative-basis states, four completed
restart states, 26 accepted electronic oligomers and six frozen fragments.
Offline replay reconstructs basis shifts and electronic decomposition; CPU
checks reconstruct 242 mode sets and both restart geometry/frequency
comparisons. Seven adverse-input checks preserve scope and identity. This
closes the named basis/restart/decomposition scientific-record gap in the
prepared package, but diagnostic artifact closure, GPU restoration, data license,
author metadata and public deposition remain release work.

The three PMDI resource pilots did not yield accepted occupancy free energies.
Two completed fixed-geometry evaluations but exceeded the unchanged 80 GiB
safety budget. The largest compact-build branch later exhausted memory in a
solvent-response product. A bounded workspace and algebraically equivalent
vector contraction pass both spin controls without changing the solvent grid,
basis, precision or physical gates. Three new same-model branches are prepared;
currently suitable H100 devices are occupied. These implementation controls do
not establish larger-system memory fitness or close M1/M4.


The subsequent immutable epoch `2026-10-23_reference_origin_review` audits the
reported DIPEA MeCN pKaH of 18.6. The author-posted primary source gives a
specific PCM-estimate note alongside a conflicting generic experimental label.
It is therefore excluded from independent measured validation. Model-to-model
sensitivity is reported separately; primary grafting sums and all-axis decisions
are unchanged. The paper and SI are rebuilt and visually checked at 14/11 pages.


A separate chemical TOC composition is now prepared from the exact PMDI
precursor and capped-product graphs, with the recorded charge ladder and an
explicitly schematic neighbor cartoon. EPS rendering, 600 dpi RGB TIFF size,
color/grayscale and small-format legibility are checked. Author review and
last-page insertion remain pending; the current paper/Overleaf export is unchanged.
The large 115-atom derivative recovery now passes four sampled total-energy
checks on a 32 GB V100; full geometry optimization is running separately. This
resource check supplies neither native G nor conformer convergence for M1/M4.

The current `2026-10-23_TOC_release_revision` epoch includes the TOC on the final
manuscript page and makes reviewer-package coverage traceable in the frozen
evidence. The checked paper is 15 pages and SI 11. The public GitHub source link
is verified; no public DOI or raw-state deposition is inferred. This supersedes
earlier notes that TOC insertion and local reviewer packaging were still pending.

## Rotor-cutoff audit and archive staging, 23 October

All 39 benchmark endpoints and 87 completed grafting reaction models now have a checked 25/50/100 cm^-1 sensitivity at fixed geometry/electronic energy/Hessian/spin. Independent full-state reconstruction checks 117 potential and 261 reaction sums. The maximum endpoint span is 25.43 mV and the maximum reaction span 4.46 kJ/mol; no nominal reaction sign changes. The revised 15-page main and 12-page SI report these settings sensitivities separately from prediction intervals and missing-conformer uncertainty. This addresses the missing benchmark rotor-sensitivity reporting; M1 and M4 remain open.

The v5 raw-state reviewer bundle is uploaded as an unpublished repository draft. Remote downloads of the archive, manifest and checksum file match local bytes. Final author/license/deposition metadata, public DOI deposition and fresh full-GPU restoration remain outstanding. The archive represents its declared earlier evidence epoch; newer cutoff/source-origin/TOC records are separately frozen in the current manuscript snapshot.

## Post-snapshot exact-PMDI numerical follow-up

The site-2 coarse force check failed at 0.001 bohr. A new component audit quantitatively attributes this discrepancy to the legacy CDS energy term. A separately predeclared full quantum fine-step test completed 15 SCFs and passed independent reconstruction without changing the 1e-6 Ha/bohr tolerance. A new electronic qualification capsule records this named validation refinement and preserves the original failure; no new optimization or active SP is claimed. Native thermochemistry and the matched PMDI cycle are still pending. This result is outside the current immutable manuscript epoch and must be explicitly described when the completed case is integrated. M1/M4 are not closed by this numerical diagnostic.

## Prior association literature, 23 October

The introduction now explicitly recognizes Wiberg et al. (J. Electroanal. Chem. 2021, DOI 10.1016/j.jelechem.2021.115224): aqueous NDI flow-battery experiments, charge-state-dependent dimer calculations, and electrolyte/steric design considerations. The checked author-repository PDF and metadata support this context citation. Its aqueous measurements do not add a MeCN benchmark endpoint or exact graft/precursor evidence. The numerical evidence epoch remains unchanged. Both documents build and the revised main/reference/TOC pages are visually checked.

## Post-snapshot response and execution follow-up

Both full spin controls pass bounded PCM first-derivative and contracted-R equivalence, and the exact 115-atom site-3 fixed-density response passes independent finite differences (maximum relative error 4.15e-6; observed GPU peak 31,224 MiB). This supports a separately gated native resource pilot, not a full Hessian/G result. These numerical controls remain outside the current manuscript numerical epoch. A subsequent shared-device audit stopped only the owned two-graft worker; its 23 evaluated checkpoints are preserved in a fresh unqualified continuation. Finite post-launch guards now inspect actual residency. M1/M4, complete native G, broader sampling and public/fresh-environment requirements remain open.

## Post-snapshot full-native allocation follow-up

Inspection showed that response-only success did not cover two later switch-curvature tensors (43.37 GiB together at the exact site-3 dimensions), so the 32 GiB native attempt was stopped rather than allowed to reach an impossible allocation. Host-backed grid tiling preserves all atoms and point ownership, and runtime bindings honor the later certified native-curvature context. Guarded full spin/Hessian replays pass to below 8e-12 atomic units. The new native pilots save an unqualified solute component and retain all physical gates, with a separately declared 31 GiB measured-peak ceiling on a minimum 32 GiB device. No PMDI native G or ensemble convergence is claimed from these controls.

## Publication narrative and fresh deployment, 23 October

The main article now leads with a chemical question and quantitative results in the abstract, introduction, discussion and conclusions. Editorial completion requests were removed from the conclusions and SI; scientific scope limits remain explicit. A new immutable deployment revision adds the separately pinned, dependency-consistent CUDA profile, full original-versus-fresh singlet/doublet controls and independent Hessian/partition reconstruction. This addresses the representative fresh-GPU reproduction requirement within its two-control scope. It does not close public deposition, cluster-independent configuration, exact-candidate occupancy or broader conformer convergence. Repository-wide consolidation remains deferred by the scientific work order.

## Spatial conformer analysis, 24 October

The four accepted doubly occupied phthalimide conformers now have exact-graph XYZ drawings, centroid/plane metrics and same-graph native G comparisons in the main paper. Both sequence spacings permit close and separated local minima. The close-contact members have lower G by 11.9/7.1 kJ/mol within these two-pose sets, including whole-chain and solvent response. Frozen coordinate data replay the drawing with original calculation paths blocked. This improves the chemical interpretation and model visualization while preserving the open M4 convergence requirement. The exact PMDI site-3 component is independently reconstructed, but the full candidate cycle remains incomplete.

## PMDI reaction-class context, 24 October

Hamilton et al. (Aust. J. Chem. 1998, DOI 10.1071/C97214) supplies primary molecular PMDI N-alkylation context. Publisher abstract and Crossref metadata are checked; full experimental procedures were not accessed. The main PMDI assessment distinguishes CsF-mediated molecular precedent, exact mono-N-ethyl precursor access and the proposed DBU/MeCN bead route. All 75 route outcomes remain unchanged; neither exact selectivity nor post-grafting redox retention is promoted. The numerical evidence epoch remains immutable.

## Current execution and exact-PMDI progression, 24 October

The public command now prefers Nitro for molecular gas-phase warm starts on supported hardware and retains the explicit V100 profile. V100 state-aware dispatch reproduces all ten frozen states/four frames; actual H100 redox parity and speed remain separate pending measurements. SI Section 8.2 documents source/cache identity and the unchanged DFT/SMD thermochemistry role. Both single-site PMDI native components reconstruct independently. The two-graft optimizer converged; its final force checks now pass after an independently controlled memory-lifecycle correction, with no repeated optimization. Canonical active energy and the matched native cycle remain pending. M1/M4 and public deposition are open.

## First-site PMDI grafting and continued conformer validation, 25 October

The new immutable first-graft evidence epoch includes both actual single-site polymer reactions and a balanced comparison against the capped local model. At 50 cm⁻¹ the separated-ion DBU reactions are −97.72/−96.96 kJ/mol, with polymer-context penalties +5.34/+6.11 kJ/mol. Common precursor/base terms cancel from those penalties. The relative site preference changes sign across the three rotor settings. This is a sampled structural comparison, not a converged ensemble or a synthesis recommendation. The 139-atom paired electronic cycle is −2.31 kJ/mol; its native G was unavailable at that first-site checkpoint. The completed neighbor-cycle section below supersedes that pending baseline status.

Broader DFT geometry refinement and canonical energy evaluation have started for the two empty-backbone configurations. Repeated final energy/force inconsistencies exclude affected conformers and have paused new geometry fan-out while a single same-energy repair pilot is tested. Completed passing geometries continue through independent energy evaluation. Scientific thresholds are retained. Runtime allocation controls and the V100-compatible tensor overlay are separately checked against full singlet/doublet Hessians; they do not establish a new polymer free energy. These execution records remain outside the manuscript's chemical narrative.

## Response qualification, stopped-search recovery and raw data, 25 October

A residual-checked GMRES response solver passes complete singlet/doublet Hessian comparisons within the original tolerance and the actual PMDI first-atom solvent-response test. A new full paired-state attempt retains every native physical gate. Its running status does not close M1 or supply a neighboring-graft free energy. Thirteen unattempted same-configuration conformers are refining under unchanged gates. Six stopped searches provide 120 exact-graph DFT warm starts; two other incomplete searches remain unprepared. These are incomplete-search seeds, not completed ensembles.

A component energy scan identifies a non-differentiable legacy-CDS carbon contact in the latest failed alternating case: one-sided slopes differ and centered curvature grows as the displacement shrinks. This changes the next scientific task from repeated optimization to validating a suitable treatment of contact minima. Failed structures stay excluded; M4 remains open. This diagnostic is retained in engineering evidence rather than added as a laboratory narrative in the main article.

The first-site polymer-context figure now replaces the composition-only capacity plot in the main paper; the latter remains in SI. The new 13-state raw-data addendum reproduces one potential and seven balanced comparisons after relocation with original checkout reads blocked. Its raw Hessians and nuclear records are included. Public deposition and general deployment portability remain open under M7.

## Reference diagnostic and solvent-surface scope

Six finished Fc single-point controls at fixed solution geometries independently reconstruct the gas/SMD electronic gaps. B3LYP and PBE0 reference-only substitutions worsen the benchmark bias; they do not provide the proposed upward reference shift or establish alternative thermochemical potentials. The article and SI report this limitation explicitly. New solvent-surface controls address a separately diagnosed contact cusp without changing primary records. Their finite molecular scope, sampled Hessian columns and unchanged physical gates are stated; no polymer G or ensemble convergence is inferred. The preceding 17-page article and 17-page SI pass frozen arithmetic, compilation and separate rendered-page review. M1, M4 and M7 remain open.

## Completed PMDI finite neighbor cycle

The 139-atom paired state and both single-site states now pass normal native admission and independent mode/nuclear reconstruction. Together with the declared shared empty reference, they give −1.07/−1.86/−1.82 kJ/mol at 25/50/100 cm^-1. The actual paired minimum has 16.30 Angstrom core separation, so the paper makes no intrinsic stacking/coupling claim. This supplies the missing baseline calculation and a source-linked PMDI example, while broader conformer ensembles remain M4’s central open requirement. The new raw control addendum passes extracted, isolated CPU replay; bulk public deposition/DOI and final metadata remain open under M7.

The [finite empty-chain subset analysis](../results/validation/empty_chain_subset_sensitivity_2026-10-26/README.md) quantifies within-set sensitivity in the alternate acetonitrile model. It does not close occupied-state ensemble convergence or enter the primary PMDI cycle.

The [proposed PMDI linker matrix](../results/discovery/pmdi_linker_matrix_2026-10-26/README.md) now binds four exact graphs to route and all-axis gates. Two already occur in the frozen library and two are new; all remain pending and no new candidate performance is claimed in the paper. Methyl shielding and cyclic-precursor access are hypotheses requiring evaluation.
