# Publication readiness: acceptance criteria and current gaps

As of 30 October 2026, **the workflow and manuscript are not submission ready**.
The objective is to address the scientific review, validate the maintained
workflow, and deliver a manuscript whose claims match its released evidence.
Code cleanup, completed example calculations and presentation updates each
contribute to that objective but do not establish completion.

The refactor moved the pre-existing manuscript, slides, calculation campaigns
and evidence to `/nfs/lambda_stor_01/homes/rzhu/archive/0_redox/`.
Its manuscript review tracker is `manuscript/review_status.md` under that root;
the machine-readable review is
`results/validation/review_comment_audit_2026-10-26/current.json`.
Those records describe the pre-refactor scientific evidence. They do not certify
the refactored checkout or establish that an archived controller is still alive.
Existing files must be preserved, and source paths/hashes must be reconciled
explicitly before a campaign is resumed.

The maintained [native-state CPU replay](../validation/publication_native_replay_2026-10-28/README.md)
now reconstructs sixty retained states and verifies 12,386 source bindings after
explicit relocation. This closes a scoped post-refactor evidence-reading gap;
it does not complete missing occupancy states or qualify current GPU execution.
Seven additional states now pass independent review and a second
[physical reconstruction](../validation/polymer_continuations_checked_2026-10-30/README.md),
verifying 10,937 source bindings. The manuscript epoch contains 67 states: 20/16
empty and 18/13 site-2 occupied states across the two backbone configurations.
The [working native collection](../validation/native_working_collection_2026-10-30/README.md)
added five more independently reviewed and reconstructed states, reaching
72: 20/18 empty and 18/16 site-2 states. Three further reconstructions now bring
the [current census](../validation/polymer_empty_batches_2026-10-30/native_census.json)
to 75, including forty empty states and a separate site-4 resource point.
The current paper integrates the component-specific empty-backbone comparison
from the 75-state frozen evidence epoch; occupied ensembles remain incomplete.
Remaining occupied graphs and complete sampling comparisons are still required.
The [remaining-graph energy roster](../validation/polymer_energy_classes_2026-10-30/README.md)
now accounts for all 160 declared inputs in the other eight graph/configuration
groups: 156 prepared energy jobs and four identical checked references. Six
geometry recoveries have explicit source relations. This closes the dependency
preparation gap, not the outstanding electronic/native calculations or ensemble
validation. A resumable resource-aware dispatcher and separate physical-review
collector carry those prepared states forward.

## Scientific acceptance

| Reviewer item | Evidence available | Requirement still to close |
| --- | --- | --- |
| M1: linked candidate decisions | Exact PMDI redox, route context, finite first-graft and four-state occupancy calculations | Carry exact enumerated designs through the linked dimensions and record evidence-based decisions; retain missing criteria explicitly. No candidate is currently demonstrated to pass every material requirement. |
| M2: redox benchmark | 39 endpoints, bootstrap intervals, grouped offset diagnostics and reference controls in the archived manuscript evidence | Check that maintained code reproduces these frozen claims. NDI accuracy and causal reference attribution remain explicit limits. The separately pinned historical calibration is not the 39-endpoint benchmark. |
| M3: base-assisted grafting | Shared per-base offset and experimental DBU/HCl and DIPEA comparisons | Keep absolute reaction signs conditional on the modeled speciation; no rank change or synthetic yield claim. Preserve the experimental-source discussion when integrating final results. |
| M4: polymer sampling | Finite primary-model cycles, two-configuration expanded starts and independently checked subsets | Complete the declared occupied-graph/configuration evaluations, duplicate review, 10-to-20 minimum/ensemble comparisons and balanced ensemble cycles. Do not mix alternate smooth-surface MeCN results with primary-model cycles. |
| M5: basis controls | Five neutral-partner geometry/Hessian controls | Retain measured potential sensitivities and state their limited scope; do not claim complete basis convergence. |
| M6: exact QC methods | Actual gradients, sampled Hessian-column checks, mode policy, spin threshold and unscaled-frequency treatment documented | Verify those checks survive the refactor and that manuscript methods describe the actual implementation. Numerical checks are distinct from external accuracy. |
| M7: public reproducibility | Cached CPU example and archived raw-state/Hessian replay packages | Reproduce the maintained workflow in a clean environment; release code, locked environments, accessible raw-state data, figure scripts, metadata and a public data DOI. Explicitly test path relocation and failure/restart behavior. |
| M8: literature | Archived audit checks 48 references, including all seven benchmark studies | Preserve the checked bibliography and verify references in the final built manuscript. |
| M9: publication narrative | Revised article/SI and a separate scientific presentation | Integrate final results into a coherent article; synchronize all tables, figures, evidence manifests, methods, limitations and data-availability statements. Build and inspect final PDFs. |

The minor conversion-constant item (m6) is addressed for the maintained pipeline:
`src/redox/core/common.py` now derives 96.48533212331002 kJ/mol/eV from the exact SI
elementary charge and Avogadro constant. The historical 96.485 factor remains
explicit for archived arithmetic; source-bound records are not rewritten.
[Migration checks](../validation/energy_conversion_2026-10-30/README.md) include
cached reaction sums and balanced conditional cycles. This units correction is
separate from experimental chemical accuracy.

## Polymer completion standard

The [two-configuration empty-backbone analysis](../validation/polymer_empty_batches_2026-10-30/README.md)
now checks all forty native input states. Chirality-preserving heavy-atom
comparisons identify 20 distinct structures for all same and 19 for alternating.
At the primary settings, expanding inputs 1–10 to 1–20 changes the finite-set
ensemble G by −1.03 and −1.23 kJ/mol, respectively. Minimum-E, minimum-G and
ensemble-G changes meet the declared 2 kJ/mol criterion across all nine
rotor/duplicate settings for both empty components. This does not complete the
occupied cycles, remove the source-search limits or prove global convergence.

Four checkpoint-recovery geometries passed independent geometry checks; their
energy and thermochemistry follow-up must use explicit source relations.
Archived native dispatch records contain completed, interrupted and upstream
pending states. Do not infer current execution from those status files: inspect
the exact host/PID/start identity and the corresponding reserved GPU first.

Local frozen-fragment electronic interactions, whole-chain conditional grafting
nonadditivity, association free energies and electronic coupling are different
quantities. The current results do not establish all four. Counterpoise controls
must pass all independent stationarity checks before their five-term sum is
reported; no gas-phase correction may silently become a solution correction.

## Maintained-workflow acceptance

1. Replay fixed, source-bound reference and candidate examples with the same
   chemical identities, charge/spin states, solvent, standard state and method.
2. Run prospective route screening before new candidate property campaigns.
3. Validate serial/parallel agreement, disjoint allocation, resume after an
   interruption, stale-owner handling and rejection of incomplete or changed
   backing records. Use measured task costs and memory to assess throughput;
   enumeration or a synthetic control is not a million-molecule QM benchmark.
4. Keep UMA gas-phase structural acceleration separate from final DFT/MeCN
   thermodynamics. Installation is not numerical qualification. The retained
   Nitro force-repeatability failures must remain visible until resolved.
5. Supply one reproducible entry point for each claimed result and figure, with
   configuration, seeds, environment, raw-record location and content hashes.
6. Rebuild the final manuscript from an immutable evidence snapshot and verify
   every quantitative claim against that snapshot. Prepare the public archive
   and DOI with the authors' metadata and data-license decision.

Completion is established by those artifacts and checks, not by a percentage
estimate or a green unit-test suite alone.
