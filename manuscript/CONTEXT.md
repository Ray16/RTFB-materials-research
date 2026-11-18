# Shared context for all redox-pipeline agents (code + manuscript)

> **Corrections (2026-11-18) — these override anything below.** The current text is
> `manuscript/manuscript_draft_v0.md`. Verified changes since this file was written:
> - Benchmark: MAE 0.247 V at the native-Hessian tier; 0.232 V at the pool's xTB tier
>   (`scripts/validation/benchmark_tier_comparison.py`). Computed with commit 46b39c3b
>   (tag `legacy/pre-refactor-main`); not re-run with the refactored code.
> - Do not justify the xTB tier with the "UKS Hessian inflates frequencies 2x" claim;
>   use the measured tier comparison instead.
> - Capacity is reported on the neutral repeat-unit mass (DMPZ ~140-170 mAh/g, anolytes
>   ~90-200); the salt-loaded values are a secondary convention (SI).
> - Pareto fronts are plain two-objective fronts with the +/-0.25 V MAE drawn, not
>   uncertainty-aware fronts.
> - lambda_i values are gas-phase (old cache schema); all DMPZ values are ranking tier.
> - QC: 269 of 271 pool states pass every ranking check; 2 flagged dianion states.
> - The scored table is verified by `scripts/validation/score_pool_legacy_tier.py`
>   (92/92 within 0.5 mV); its pipeline audit stage had failed.
> - Bias-corrected proposed cells: about 1.7-1.9 V. Tests: 780 passing.


Repo root: /nfs/lambda_stor_01/homes/rzhu/0_redox_refactor
Archive (prior manuscript, slides, validation campaigns, bulk calcs — READ ONLY,
do NOT copy back into the repo): /nfs/lambda_stor_01/homes/rzhu/archive/0_redox
Target journal: **JCIM (Journal of Chemical Information and Modeling)**.

## Hard constraints (from AGENTS.md + CLAUDE.md) — every agent obeys
- NO GPU work. Do NOT launch DFT, Hessian, UMA, or any QM. Do NOT set
  CUDA_VISIBLE_DEVICES. A campaign (`results/expand_catholyte_20261117`) is LIVE
  and the cluster is shared. CPU-only validation (pytest) is fine with
  OMP_NUM_THREADS=1 / MKL_NUM_THREADS=1 / OPENBLAS_NUM_THREADS=1.
- Do NOT touch certificate-locked files: src/redox/qm/qrrho.py,
  src/redox/qm/xtb_kernel.py, src/redox/qm/xtb_runtime.py,
  src/redox/quality/sampling_snapshot.py, and the hash-locked example payloads.
  Reformatting them invalidates the projected-thermal certificate.
- Keep half-reactions and reorganization components as SEPARATE quantities.
  **No weighted-sum score.** Use real chemical names. A computed descriptor or a
  completed job does NOT establish synthetic feasibility or material performance.
- A capped graft, a soluble parent, and an oligomer are DIFFERENT models.

## The pipeline (scientific spine) — manuscript Methods source of truth
Per-candidate driver: `redox.workflows.motifs.commands()`. Nine stages:
  build → UMA pre-opt → DFT geometry+SMD(MeCN) per state [neu/ox1/ox2]
  → single-points → thermochemistry → properties → QC coverage/cross-points
  → scoring/Pareto.
- Level of theory (CORRECTED from code — composite protocol): geometry optimized
  at **r²SCAN-D4 / def2-SVP(D), SMD(MeCN)**; energy scored at **ωB97M-V(VV10) /
  def2-TZVPD, SMD(MeCN)**, RI-J, conv_tol 1e-9. (ωB97X-D3 was planned but is not
  implemented in gpu4pyscf → ωB97M-V substituted.) ε_s=37.5, ε_op=1.806, n=1.344;
  298.15 K; potential window [-3.0, 2.0] V vs Fc, anolyte/catholyte divider 0 V.
- Grafting model: SMILES motif + 4-methylbenzyl (Merrifield proxy) cap; redox
  states {neutral, ox1 (+1), ox2 (+2)} (or reductions for n-type anolytes).
- Thermochemistry is TWO-TIER. CORRECTIONS from code:
  * The honest rationale for the cheap uniform RRHO tier is NOT "basis
    consistency" (xTB-RRHO has no basis). It is: (a) RRHO corrections are
    ~method-transferable so a uniform cheap Hessian cancels in energy DIFFERENCES,
    and (b) gpu4pyscf's UKS analytic Hessian is broken (inflates open-shell freqs
    ~2×), making DFT Hessians unsafe for the radical states. USE THIS RATIONALE.
  * "native DFT qRRHO" is loose: native.py uses a DFT (r²SCAN-D4/def2-SVP(D)/SMD)
    Hessian but evaluates the qRRHO partition with the certificate-locked PROJECTED
    internal-mode xTB 6.7.1 kernel (tag native_solution_dft_qrrho_spin).
  * The existing 91-candidate pool CSV used the LEGACY xTB-RRHO tier
    (full_pool_scored_legacyG_...). But the current pipeline CODE default is
    thermal_definition="native" for fresh runs, and native thermochem is GPU-only
    (refuses CPU). So "legacy = pool-wide ranking" describes the existing table,
    not the current default. Frame as deliberate two-tier design.
- Reference: potentials vs ferrocene (Fc). Electrolyte: config encodes counterion
  PF6- + window + divider; TBAPF6 is the implied supporting salt in MeCN.

## Reproducibility assets (Methods / Data-availability)
- 5 byte-exact portable replays (examples/), 744 passing contract tests,
  certificate-locked qRRHO kernel. See README.md, MIGRATION.md, docs/.

## Key quantitative results (VERIFY against results/ + archive before final use;
## mark anything you cannot confirm as [TODO: verify <source>])
- Pool: **91 candidates** scored
  (results/validation/full_pool_scored_legacyG_MeCN_TBAPF6.csv).
- Composition: **85 anolyte / 6 catholyte**; 89/91 are 2-electron couples.
- Families: imide (n-type) 36, quinone (n-type) 29, pyridine-multi-e 18,
  phenazine (p-type) 5, mono-imide (n-type) 2, phenothiazine (p-type) 1.
- Validation benchmark (RECONCILED — use these, the brief's "0.226 V / 44
  couples" was WRONG): **MAE 0.247 V, RMSE 0.279 V, bias +0.19 V, max|err| 0.547 V,
  over 39 tier-A endpoints (redox couples) from 30 molecules and 7 studies**;
  macro-MAE 0.263 V. Source of truth: archive
  manuscript/evidence/2026-10-22_benchmark_complete/benchmark.json, corroborated
  in archive manuscript.tex and docs/PUBLICATION_READINESS.md. (One E[event] =
  one redox wave = one couple = one endpoint; 39 tier-A admitted of 47 curated.)
  The "0.226" in the brief was a coincidence — it is one endpoint's
  serialized_calc_V, not the MAE. Wave-1 bias +0.266 V, wave-2 bias -0.064 V →
  bias is systematic → RANKING more reliable than absolutes.
  M2 CAVEAT STILL OPEN: these are frozen-archive values; re-confirming the
  MAINTAINED/refactored code reproduces them requires a QM re-run (deferred under
  the no-GPU constraint). Do not state "validated by the maintained code" until
  that re-run is done.
  !! The Introduction (01) and Abstract/Conclusions (00) drafts were written with
  the WRONG 0.226 V / 44 couples — fix to 0.247 V / 39 endpoints during stitch.
- Catholyte gap / DMPZ (5,10-dihydrophenazine) 2e- catholytes, E vs Fc (from
  the scored CSV, VERIFY): N-mb-N'-ethyl +0.43; N'-methyl +0.47;
  N'-(2-methoxyethyl) +0.47; N'-(2,2,2-trifluoroethyl) +0.62;
  2,7-dimethyl-N'-ethyl +0.30; and a 3,7-dimethoxyphenothiazine +0.78.
  Design rule: electron-withdrawing N'-substituent raises E, ring-methyl lowers E.
- Proposed full cells ~1.9 V (DMPZ catholyte x quinone/imide anolyte) — PROPOSAL
  only; cell voltage = E(catholyte) - E(anolyte). Not demonstrated capacity/
  cycling/crossover/synthesis.

## Property-confidence disposition
- Potentials vs Fc: strong, validated → main text.
- Usable capacity (TBAPF6/MeCN): strong → main text (Pareto axis).
- λ_o (outer-sphere): defensible continuum estimate (Born + molecular-cavity PCM).
- λ_i (inner-sphere): CURRENTLY low-confidence (gas energies on SMD geometries,
  single-conformer). A fix is in progress (all-SMD self-consistent surface +
  matched-conformer averaging). Until that lands, treat λ_i as SI/qualitative.
- Scoring: σ-aware Pareto fronts, anolyte & catholyte pools SEPARATE. No scalar score.

## Shared planning docs (read these)
- manuscript_outline.md — figure-by-figure spine (method-led).
- refactor_plan.md — execution/parallelization plan (scheduler, tables, ferrocene).

## Coordination
- Code agents each work in an ISOLATED git worktree on disjoint files; produce a
  committed branch + a diff summary; run the relevant CPU tests. No GPU.
- Manuscript agents each write ONE section file to
  scratchpad/manuscript/<section>.md. Do not edit another section's file.
  JCIM house style: concise, third person, SI for heavy tables/derivations.
