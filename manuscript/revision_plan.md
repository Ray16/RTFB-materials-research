# Colleague review: revision ledger

Review received 22 October 2026, referring to the immutable 21 October draft.
The review's scientific concerns remain applicable even where counts have since
changed. A task is closed by checked evidence and a corresponding manuscript
revision, not by a proposed calculation or a running job.

## Major comments and submission requirements

| Comment | Current evidence and required revision | Status |
| --- | --- | --- |
| M1: no end-to-end decision | Exact enumerated PMDI and NDI matches now have freshly checked six-state redox/Fc reconstruction, grafting hashes and two mapped route assessments (both hold). This is a partial join; complete missing route and occupancy assessments. Present a decision table showing which evidence changes priority or leaves synthesis on hold. Method controls must not manufacture synthesis eligibility. | Open; highest scientific priority |
| M2: benchmark bias/reference | Keep all 39 unfitted endpoints primary. Report the existing molecule-cluster intervals, qualify the seven unequal study clusters, and add endpoint-, molecule- and study-held-out offset analyses, including redox-step-specific models. Compare Fc-independent spacings and multiple organic-reference diagnostics. Verify the Fc literature's methods and conventions. Discuss the orbital/Jahn–Teller limitation without assigning an unvalidated correction. | Analysis checked and included; reference attribution remains conditional |
| M3: constant base shift | State explicitly that a fixed base representation adds the same constant to each reference reaction and cannot change within-base rankings. An exact-graph MeCN table check now gives DBU pKaH 24.31 and HCl 10.30, a 22.81 kJ/mol discrepancy for the capture term. This does not uniquely locate the error because HCl cancels in the full cycle. DIPEA remains conditional; do not transfer DBU's discrepancy as a DIPEA correction. Four additional component calculations test two HCl/Cl-free relative-basicity reactions. | Two relative controls independently checked and included; source/speciation qualification open |
| M4: conformer/tacticity sensitivity | All six selected phthalimide G values now pass, giving -5.15/-11.23 kJ/mol, but this does not establish sampling convergence. Document seed 20261020, 4 embeds, 2 xTB seeds and 2 DFT poses per graph. Add broader conformer search, distinct backbone stereochemical scenarios, DFT re-ranking and a defined ensemble treatment. Predeclare convergence criteria. Until demonstrated, retain the result as a finite method-control example, not a spacing design rule. | Open; highest scientific priority |
| M5: geometry/Hessian basis imbalance | Five matched neutral-SVPD controls independently reconstruct, giving -11.7 to -35.1 mV shifts with the existing anion/Fc fixed. The main paper and SI now report geometry, nuclear-partition and resulting potential changes separately. Keep primary source-frozen records unchanged. Alternative-basis controls must pass their physical checks and cannot inherit the original certificate by relabeling. | Five-control sensitivity checked and reported; full basis convergence not claimed |
| M6: Methods/QC precision | State exactly two tested Hessian columns and both displacement sizes; explain |nu| for tolerated -5..0 cm^-1 modes, unscaled frequencies, spin-contamination threshold 0.1, actual force values and xTB symmetry detection. Remove the unsupported adjective “conservative” for the 4.5e-4 force threshold. Cite qRRHO and program/method papers. | Exact gate description included; further physics controls remain separate |
| M7: portable/public reproduction | A fresh isolated CPU environment now reproduces graphs, route holds and a figure; that earlier example did not reproduce a potential or grafting G from raw states. A relocated standard-library example now reconstructs a potential and both HCl-reference/DBU grafting sums from ten hashed raw states. Explicit conda artifacts and pip versions are recorded; xTB/build and core scientific dependencies are pinned. Portable primary paths, fresh GPU installation and a full raw-state/public DOI release remain open. Prepare a reviewable public archive and deposition metadata; a DOI/public deposition is not yet created. Retain strict integrity checks for primary scientific admission. | Partly implemented; full release open |
| M8: literature context | Verify each DOI and source. Add computational electrochemistry/reference, high-throughput design, redox-polymer, workflow/provenance and missing method literature where it supports a specific argument. Do not turn remembered references into unverified citations. | Expanded from 14 to 27 verified cited references; further relevant context in progress |
| M9: cohesion and hedging | Lead with the chemical decision enabled by the identity-linked framework. Move synthetic scheduling, concurrency, restart and numerical-repair detail to a separately built SI. Consolidate limitations, define unfamiliar terms and reduce implementation vocabulary in the main article. | Main and separate SI built; 12 and 10 pages, respectively |
| Remaining endpoint | The fixed panel now has 39/39 endpoints; the last triplet alternative is 0.707 eV above the selected singlet electronically. Independent reconstruction: 111 checks, zero discrepancies. Alternative-spin thermochemistry remains explicitly unresolved. | Calculation/reconstruction resolved; latest draft updated |

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
| m8: figures/graphical abstract | Replace the box-only abstract with a chemically informative composition and move implementation/status banners to captions. Check current JCIM requirements. |
| m9: jargon | Replace implementation labels where possible and define the remaining evidence/decision terms in one compact table. |
| m10: reproducibility details | Record seed, sampling counts, backbone configuration and xTB symmetry detection explicitly. |
| Hessian replay rounding | Correct 1.34e-11 to the value supported by the stored evidence (1.333e-11; 1.33e-11 at three significant figures). |
| QC/imaginary-mode treatment | Match the text to actual sampled columns, thresholds and tolerated-negative-frequency handling. |
| Base-induced sign changes | Show the per-base shift across the whole qualified panel, including AQ and NDI, without a ranking-change claim. |
| xTB build/path locking | Pin and document the actual quantum environment; configure deployment paths without silently bypassing scientific integrity or admission. |

## Work already running when this review arrived

The stacked-trimer optimizer returned convergence and an active SP, then a later library call crashed. Separate fixed/fine final-force checks failed, so the case remains excluded. Its obsolete finite-panel monitors have stopped; the two broader CREST searches continue independently. Its small-step force diagnostic is complete. They remain distinct from the completed six-state spacing
panel. The two-spin quantum interruption experiment has completed and passed an
independent reconstruction of its four final thermal partitions. The four relative-basicity components and both balanced transfer sums have passed independent reconstruction. The two original empty-backbone pilots remain live. Twelve full-parallel searches now cover all six occupancy graphs at two configurations, including two reference repeats; actual starts and disjoint CPU sets are verified. DFT re-ranking and ensemble G remain pending.
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

## Conformer thermochemistry follow-up, 23 October

A finite partition implementation and an all-pose native thermal adapter now
reconstruct the six inherited spacing G records against their active electronic
and projected nuclear/spin parents. The other six original poses are requested
in a new source-bound campaign; a reserved 119-atom native pilot has actually
started. This quantifies the two-pose baseline while wider CREST searches run.
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
