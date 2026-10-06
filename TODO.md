# TODO

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
- [x] Figures `results/figures/candidates/discovery_*.png` (capacity vs SA, chemical-space PCA,
      D3TaLES lambda; only physical D3TaLES lambda 0-1.5 eV is used — 58 of 113).
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
- [ ] `plot_lambda_decomposition.py` still reads the preliminary Born λ_o
      (`results/lambda_outer_merrifield_born.json`) while `results/lambda_outer.csv` (PCM v2) is
      canonical — switch it to `lambda_outer.csv`.
- [ ] `scripts/polaris/polaris_build_redox_env.sh` points at the expired `/eagle/FoundEpidem`
      allocation — repoint to `/grand/FRAME-IDP` before the next Polaris run.
