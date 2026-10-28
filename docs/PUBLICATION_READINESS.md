# Publication readiness: acceptance criteria and current gaps

As of 28 October 2026, **the workflow and manuscript are not submission ready**.
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

The remaining minor conversion-constant item (m6) also needs a controlled
migration: `src/redox/core/common.py` still uses the historical 96.485 kJ/mol/eV.
New examples use the exact SI factor and label the historical comparison.
Changing a constant must not silently invalidate existing source-bound evidence.

## Polymer completion standard

The twenty-conformer empty-backbone analysis is at
`results/validation/empty_chain_twenty_2026-10-27/` in the archive. For one
configuration, expanding ranks 1–10 to 1–20 leaves the minimum E and G unchanged
and changes the finite-set ensemble G by −1.03 kJ/mol. This meets the declared
2 kJ/mol batch diagnostic **for that component only**. It does not complete the
occupied cycles or prove global conformer convergence.

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
