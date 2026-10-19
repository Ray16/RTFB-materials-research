# TODO

## Current priorities (2026-10-19)

- [x] Integrate source-bound two-spin equilibrium in both readers and ranking gates;
  transactional selection verified on CN-PMDI, with degeneracy counted once.
- [x] Extend integrity to every selected spin component and independently rebuild
  bond graphs. Keep mixed-spin kinetic descriptors withheld; 428 tests pass.
- [x] Rebuild and check equilibrium/candidate exports and the 39-page slide PDF;
  436 software tests, 30,455 global checks and 572 candidate checks pass.


Use [LATEST_REVIEW.json](results/LATEST_REVIEW.json) and the
[candidate feed](results/candidates/current.json) for current counts. Sections
below this checkpoint preserve historical worklists and superseded estimates.

- [x] Integrate independently qualified internal-mode thermochemistry, including
  genuine soft modes, into the producer and both primary readers; correct matched Fc.
- [x] Correct the geometry driver's grid response and optimizer-convergence reporting;
  391 software tests and an actual CPU geometry optimization pass.
- [x] Make generic native repairs carry corrected partitions through validated
  atomic promotion with rollback; preserve convergence on resumed optimizations.
  The full suite now passes 395 tests.
- [x] Independently review and promote the NDI methoxyethyl and NDI butyl radical
  minima (157 checks each, zero failures), preserving exact original-state archives.
- [x] Complete corrected potential pairs for 24 of 25 registered new candidates;
  572 export arithmetic checks pass with zero discrepancies in the current snapshot.
- [x] Finish a stable corrected global property-table review, refresh slides and
  synchronize the result index with verified snapshots.
- [~] Finish hexyl NDI's native surface under unchanged physical gates. Both prior
  coarse-Hessian failures are preserved; one qualified finer-step H100 replay is
  verified running and the second input waits for suitable free capacity. CN-PMDI's
  equilibrium and ring-methyl NDI's reviewed minimum are now admitted.
- [ ] Broaden matched experimental anchors across viologen, PMDI/NDI, quinone and
  positive-electrode families. Track solvent/reference/structure identity and
  quantify potential, wave-spacing and conformer uncertainty. Seven potentials
  from four parent molecules do not establish exact-graft accuracy.
- [ ] Qualify solution reorganization and kinetic descriptors with controlled solvent
  response and conformer matching; keep gas relaxation diagnostics separate.
- [ ] Assess synthesis, tether stability, charge-state-dependent chemical degradation,
  solubility/swelling and transport for the shortlisted exact structures. Choose
  targeted experimental tests where computation cannot establish performance.
- [ ] Continue balanced candidate exploration using complementary curated sources,
  exact-graph deduplication and matched substitutions. Use uncertainty and chemical
  diversity to select the next batch rather than only nominal voltage/capacity.
- [ ] Finish reproducibility: document a small end-to-end example, record the
  validated environment, and publish one source-consistent set of tables, figures,
  slides and field-specific QC. Keep reusable workflows in the package and dated
  scientific evidence separate from supported entry points.

## Validation and discovery continuation (2026-10-10)

2026-10-11 continuation:
- [x] Recover and verify Fc+ gas ground branch; production occupation seeds prevent repeat error.
- [x] Replace `d3_80izcq/red2` triplet with1.757eV-lower singlet; archive every old file, retain triplet helper and verified spin override.
- [x] Complete3 nonequilibrium solution controls; preserve clear total/inner conventions and diagnostic status.
- [x] Complete AQ local/linker-fixed controls; validate GFN2-ALPB minima by mode following for2 targeted states.
- [x] Correct mixed gas/SMD reference diagnostics and charged-as-neutral solvation fallback.
- [x] Gate unstable RRHO standard properties; retain estimates separately. New checkpoint:10,223 checks, zero discrepancies.
- [x] Register2 new motifs and verify5 converged, chemically intact UMA warm-start states.
- [~] Full conservative-gradient solution Hessian; integration-grid response defect diagnosed, numerical reference running.
- [~] Complete all remaining UMA<0.5eV DFT spin checks; electronic margins for3 original singlets verified.
- [ ] Validate transferable thermal protocol, restore corrected d3_80izcq G/cross records and regenerate final tables.
- [x] Preserve original results; regenerate a reviewed snapshot in `results/review_2026-10-10_final/`.
- [x] Independent audit expanded to multi-electron summary and workbook: 8,306 checks, zero discrepancies.
- [x] Exclude partial kinetics from complete Pareto objectives; require consecutive reductions from resting state.
- [x] Validate new motif attachments, stereo, charge/spin parity and ladders; prevent ID collisions and unsafe rebuilds.
- [x] Automate motif-to-audited-table execution with reservations, writer locks, CPU caps/affinity and stage logs.
- [x] 48 contract tests, real new-motif embedding/resume test, cached viologen integration run; environment verified.
- [x] Check Wang SI Table S1: demote ambiguous fluoro-benzoquinone E2/spacing to tier B.
- [x] Refresh discovery without duplicate grafted molecules: 104 new quinones, ten preliminary-front molecules.
- [x] Generate nine uncomputed PMDI/NDI proposals with validated grafting and state definitions.
- [ ] Repair 11 incomplete declared states; retain the 15 imaginary-mode thermal review flags until resolved.
- [ ] Compute solution/conformer-matched reorganization for flagged steps; DFT spin checks for small UMA gaps.
- [ ] Review grafting chemistry of the 19 next proposals, then register selected motifs and run production QM.

Details, scientific limits, commands and artifact links: `docs/VALIDATION_REVIEW_2026-10-10.md`.
Parallel continuation: `docs/PARALLEL_PROGRESS_2026-10-10.md`.
- [x] Recover ten missing SP records by proven adoption; source geometries/results unchanged.
- [~] Reserved lambda1 repair: one missing DFT state and 28 cross points; bounded finalizer will audit a new snapshot.
- [x] Chemical identity, unit and SCF-selection audit; flag decomposed `mv_ip1/s0` helper INVALID.
- [x] Fix full cross-cache provenance, all-leg anion QC, radical/torsion mapping and thermal cache/reference contracts; 58 tests pass.
- [x] Parallel curated expansion: 39 proposals (38 connectivities), 156 independent descriptor checks pass, all structures embed and build-check.
- [ ] Close remaining thermal/spin/solution-lambda physical validation limits before claiming every metric is validated.

Living worklist (strategic plan: `docs/PLAN.md`; axis status: `docs/DESIGN_AXES.md`).
Legend: `[ ]` todo · `[~]` in progress · `[x]` done · `[!]` blocked. Previous list:
`archive/docs/TODO_2026-08.md`.

## Data contracts (review round 2, 2026-09-29)
- [x] Integrity gate: missing verdict = INCOMPLETE, not REJECTED (had rejected 4 candidates).
- [x] Rename `reversible` -> `intact_bound` (`redox.properties.integrity`, `results/state_integrity.csv`).
- [x] Fc/Fc+ reference: same G definition as molecules, computed LIVE at the active protocol.
- [x] Uniform energy basis: `redox.core.protocol.ACTIVE_SP` (def2-TZVPD for all states + Fc);
      protocol-addressed energy records (`calcs/dft/<id>/<state>/sp/<hash>.json`).
- [x] Reorg cross-point caches protocol-addressed; `_cache_ok` checks xc/basis/NLC/convergence.
- [x] `dft.run_batch` flags a result.json from another optimization protocol as STALE.
- [x] Thermal gate: missing thermal -> INCOMPLETE (no silent 0); imag modes flagged.
- [x] Explicit lambda conventions (lambda_O, lambda_R, sum, lambda_het, lambda_se); removed
      the mixed `lambda_total`.
- [x] Pareto: complete rows only on the primary front; per-row sigma columns.
- [x] Contiguous per-pool capacity from the resting state; one scorecard row per pool.
- [x] lambda_o: sign bug masked by abs() fixed; point-charge sphere test == Born; SCF check;
      diffuse default; metadata; both geometries.
- [x] Fail-fast finalize chain incl. stability + lambda_o + tests; pytest suite `tests/`.
- [x] Uniform-basis SP campaign (pipeline states + cross points + OROP): 173 + 348 records, all
      with proven provenance (2026-10-02).
- [x] lambda_o v2 for all 21 candidate couples — now actually consolidated into `results/lambda_outer.csv`.
- [~] Thermal-model sensitivity (sthr 25/50/100) -> `results/thermal_sensitivity.csv`.
- [~] DFT spin check aq_benzylamino dianion triplet (CPU single point).

## Accuracy of computable axes (before screening)
- [~] E° offset: level crossing at OROP's own level (B3LYP-D3/6-31G*, C-PCM) shows our organic
      absolute potentials agree with OROP's (~0.07 V) and the gap sits in the Fc reference
      (ours 4.185 V vs OROP's tabulated 4.662 V) — FINDINGS #24. Remaining: why OROP's Fc number
      is higher (SI Text S2 not accessed); the family-dependent E1 offset on the sourced
      benchmark (quinones +0.31, imides +0.12 V). Until explained: rank within family only.
- [ ] Per-family E° calibration once (above) is resolved; needs >=3 anchors per family
      (imide family has 0 experimental anchors — add PMDI/NDI literature E°).
- [x] Second-reduction accuracy: 8 sourced tier-A E2 (benzoquinones, 1,2-NQ, PMDIs): bias −0.07 V,
      SD 0.18 V. Open: 1,2-naphthoquinone −0.46 V outlier; no anthraquinone E1/2 second wave yet.
- [x] Disproportionation benchmark n=3 -> 10 sourced tier-A spacings: ΔG_disp over-estimated by
      ~24 kJ/mol (FINDINGS #24). Next: explain the E1 offset (E2 is accurate).
- [ ] Dianion lambda: gas-phase 4-point is ill-defined when the gas dianion is unbound
      (HOMO>0) — score the 2nd reduction's lambda in SMD or flag as not computable.
- [ ] Spin: DFT-check every state whose UMA gap < 0.5 eV (UMA underestimated 0.19 vs 0.81 eV).
- [ ] lambda_o: validate PCM vs a trusted nonequilibrium implementation (Q-Chem/Gaussian) on
      2-3 molecules; the point-charge sphere test only validates the operator path.
- [ ] Conformer-matched core lambda for `conformer_jump` couples (tether rotation, not core).
- [ ] Stability beyond disproportionation: per-family decomposition library (isodesmic dG).
- [ ] Commit + push this round once the campaign lands (user go-ahead needed).

## Validation round 3 (2026-10-02) — FINDINGS #24
- [x] Independent recompute audit of every published table (step 12 of finalize): 0 mismatches.
- [x] SCF ground-state check of every candidate energy (up to 9 guesses); open shells now keep
      the lowest of minao/atom/huckel.
- [x] Provenance: records adopted from pre-RI-J result.json recomputed (7-12 meV); adoption
      now requires recorded density-fitting / tolerance / guess-sweep provenance.
- [x] Sourced MeCN E1/E2/spacing benchmark (42 literature rows, 33 new molecules computed).
- [x] sigma_E / sigma_disp from grounded points only; Fc no longer counted as a validation point.
- [ ] Explain the E1 offset (quinones +0.31 V, imides +0.12 V; E2 accurate) — continuum ion
      pairing of the radical anion? Fc reference (gas IE 0.29 eV low vs NIST)? Do not fit it.
- [ ] Ferrocenium gas-phase SCF (2.5 eV high; only an occupation guess finds the minimum):
      recompute if any Fc gas quantity is ever needed.
- [ ] Verify Pavlishchuk & Addison 2000 constants from the paper (SCE/Ag+ -> Fc); until then
      values converted with them stay tier B.
- [ ] Anchors without a verified source: phenothiazine (+0.26), anthraquinone (-1.28/-1.90),
      N-methylpyridinium (-1.8) — `grounded=False`; find primary MeCN data or drop.

## Candidate discovery (2026-10-06)
- [x] `redox.screening.discovery` + `scripts/mining/identify_candidates.py`: rules regenerate
      16/16 registered candidates; D3TaLES -> 113 new graftable 2e- quinones (107 p, 6 o);
      0 graftable bis-imides (3 found, all N,N'-dialkyl), 0 viologens (none in D3TaLES).
- [x] Figure `results/figures/candidates/discovery_capacity_vs_sa.png` (pool vs registered
      candidates on identical axes); only physical D3TaLES lambda (0-1.5 eV) is used.
- [ ] Stage 2: run the 8-molecule pre-filter front (`results/discovery/d3tales_prefilter_front.csv`)
      through the full grafted pipeline; check halo-quinones (Cl/Br on the quinone ring) for
      nucleophilic substitution — the grafting conditions themselves may attack them.
- [ ] Design generators for bis-imides and viologens (not available from D3TaLES).
- [ ] Validate a cheap stage-1 ranking (UMA/xTB) against our DFT before cutting deeper.

## Repo hygiene (refactor 2026-10-02)
- [x] `src/redox` split into subpackages `core/ build/ qm/ properties/ screening/ validation/`
      (`redox.redox` -> `redox.properties.potentials`); `scripts/` regrouped (pipeline, analysis,
      validation/{orop,solvation,reorg_d3tales,reorg_literature}, mining, plotting/<figure dir>);
      obsolete launchers/GPU enforcers/dead modules -> `archive/` (see `archive/README.md`).
- [x] `run_uma.sh` / `run_dft.sh` reserve GPUs via `gpu_reserve run` (no hand-pinned CUDA ids).
- [x] `plot_lambda_decomposition.py` (preliminary Born λ_o) archived 2026-10-06 with the other
      one-off figures (archive/README.md); λ_o in the scorecard is the PCM v2 table.
- [ ] `scripts/polaris/polaris_build_redox_env.sh` points at the expired `/eagle/FoundEpidem`
      allocation — repoint to `/grand/FRAME-IDP` before the next Polaris run.
