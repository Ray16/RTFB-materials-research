# Manuscript outline — method-led spine

**Working title:** *A validated, reproducible pipeline for screening graftable
redox motifs for non-aqueous organic redox flow batteries*

**One-sentence thesis:** We present a validated, end-to-end computational
pipeline that builds exact grafted-motif models, computes reliability-gated
redox properties referenced to Fc, and ranks candidates on a transparent
multi-property front — and we demonstrate it on 91 candidates, where it
surfaces grafted dihydrophenazine 2e⁻ catholytes to fill an otherwise
anolyte-dominated design space.

**Spine (method is the contribution; molecules are the demonstration):**
problem (screening graftable motifs needs *exact* models + reliability gating +
external validation, which ad hoc DFT screens lack) → the 9-stage pipeline →
validation against experiment (the credibility anchor) → applied to 91
candidates → showcase result (the catholyte gap, filled with DFT-backed design
rules) → proposed full cells.

Contract with AGENTS.md: keep half-reactions and reorganization components as
SEPARATE quantities; **no weighted-sum score**; real chemical names; a computed
descriptor ≠ synthetic feasibility or material performance.

---

## Figure / section plan

Each item: **Claim → matched comparison → evidence (status) → figure.**

### Fig 1 — The grafting abstraction (schematic). *Methods, framing.*
- **Claim:** redox motifs grafted to a preformed Merrifield-type backbone can be
  modeled exactly by a 4-methylbenzyl (Merrifield proxy) cap, with redox states
  defined per motif.
- **Evidence:** `structures/grafting.py`, `structures/attachment_identity.py`;
  912-graph identity registry. STATUS: solid.
- **Figure:** schematic — SMILES motif + cap → 3D → {neu, ox1, ox2} states.
- Caveat line to include: capped graft ≠ soluble parent ≠ oligomer (three
  distinct models), stated up front (README already says this).

### Fig 2 — The pipeline (flow diagram + method table). *Methods, the contribution.*
- **Claim:** a single reproducible pipeline takes a motif to reliability-gated
  redox properties: build → UMA pre-opt → ωB97M-V/def2-TZVPD + SMD(MeCN)
  geometry per state → single-points → thermochemistry → properties → QC gates.
- **Evidence:** `redox.workflows.motifs.commands()`; level of theory in
  `backends/dft.py`, `backends/sp.py`; byte-exact reproducibility (5 replays,
  744 contract tests, certificate-locked qRRHO kernel). STATUS: solid.
- **Figure:** 9-box flow diagram with the method table (level of theory,
  solvent, standard state, reference) as an inset; annotate which stages are
  GPU/DFT vs CPU.
- **Method caveat (own subsection, not hidden):** thermochemistry has two tiers
  — legacy xTB-RRHO (pool-wide, basis-consistent → used for *ranking*) vs native
  DFT qRRHO (disk-gated, expensive → used to *confirm leads*). Frame as a
  deliberate two-tier design, not a weakness.

### Fig 3 — Reliability gating. *Methods: why the numbers are trustworthy.*
- **Claim:** every reported candidate passes explicit QC: spin/⟨S²⟩, geometry
  integrity, full state coverage, and cross-point consistency; candidates that
  fail are held, not silently dropped.
- **Evidence:** `quality/readiness.py` (coverage, --strict), `quality/spins.py`,
  `quality/stability.py`; the `status`/`coverage` stages. STATUS: solid.
- **Figure:** QC gate funnel (N registered → N state-complete → N QC-clean →
  N scored) OR a small table of gate pass-rates. Reinforces the method claim.

### Fig 4 — Validation against experiment (THE CREDIBILITY ANCHOR). *Core result 1.*
- **Claim:** pipeline potentials reproduce experimental redox potentials with
  **MAE 0.226 V and a +0.19 V systematic bias** across **44 couples**.
- **Matched comparison:** computed E vs Fc vs literature E vs Fc, same solvent,
  same reference, bootstrap intervals, grouped-offset diagnostics.
- **Evidence:** `quality/benchmark_panel.py`, `reporting/benchmark.py`,
  `properties/potentials.py`, pinned Fc reference. STATUS: solid (per prior
  archive benchmark; REQUIREMENT: re-confirm the maintained code reproduces the
  frozen 39-endpoint/44-couple numbers — this is reviewer item M2).
- **Figure:** parity plot (computed vs exptl, colored by family) with MAE/bias
  annotation + residual histogram. This figure is what buys the reader's trust
  in Figs 5–6.
- Honesty line: the +0.19 V bias is systematic → ranking is more reliable than
  absolutes; NDI/causal attribution stated as explicit limits.

### Fig 5 — Applied to 91 candidates: the multi-property front. *Core result 2.*
- **Claim:** across 91 candidates the pipeline maps the accessible design space;
  the space is strongly **anolyte-dominated (85 anolyte / 6 catholyte)**, with
  families: imide (36), quinone (29), pyridine-multi-e (18), phenazine (5),
  mono-imide (2), phenothiazine (1).
- **Matched comparison:** σ-aware Pareto fronts, anolyte and catholyte pools
  separately (potential × usable capacity in TBAPF₆/MeCN), NOT a scalar score.
- **Evidence:** `screening/pareto.py`, `screening/scorecard.py`,
  `properties/capacity_and_proxies.py`; `full_pool_scored_*.csv`. STATUS: done.
- **Figure:** two Pareto panels (anolyte | catholyte), points = candidates,
  error bars from σ, families by color/marker. The visual sparsity of the
  catholyte panel sets up Fig 6.

### Fig 6 — SHOWCASE: filling the catholyte gap. *Core result 3 (the payoff).*
- **Claim:** grafted **dihydrophenazine (DMPZ)** 2e⁻ catholytes populate the
  empty catholyte front; substituent tuning shifts E vs Fc predictably
  (e.g. ethyl ≈ +0.43 V → methoxyethyl +0.47 → trifluoroethyl +0.62;
  2,7-dimethyl pushes down to +0.30; methoxy/phenothiazine up to +0.78).
- **Matched comparison:** matched DMPZ series differing by ONE substituent
  (the `expand`/`expand_ring_panel` matched-hypothesis machinery); read the
  design rule off the matched deltas.
- **Evidence:** `screening/expansion.py`, `config/custom_candidates.py`, active
  `results/expand_catholyte_20261117` run. STATUS: in progress (4 DMPZ
  catholytes computing now) → this figure is the one still being filled in.
- **Figure:** DMPZ scaffold with substituent positions annotated by their
  E-shift (design-rule figure) + the catholyte Pareto panel re-shown with the
  new points landing in the gap.

### Fig 7 — Proposed full cells. *Implication / outlook.*
- **Claim:** pairing the surfaced DMPZ catholytes (~+0.4–0.6 V) with pipeline
  anolytes (quinone/imide, ~−1.0 to −1.5 V) gives proposed non-aqueous full
  cells near **~1.9 V**.
- **Matched comparison:** cell voltage = E(catholyte) − E(anolyte), both from
  the same validated pipeline/reference.
- **Evidence:** derived from Figs 5–6. STATUS: proposal.
- **Figure:** cell-voltage ladder / pairing diagram.
- **Hard caveat (required):** a computed cell voltage is a proposal, NOT
  demonstrated capacity, cycling, crossover, or synthetic feasibility.

---

## Property-confidence disposition (what goes where)
- **Potentials vs Fc** — strong, validated → main text (Figs 4–7).
- **Usable capacity (TBAPF₆/MeCN)** — strong → main text (Fig 5 axis).
- **λₒ (outer-sphere)** — defensible continuum estimate (Born + molecular-cavity
  PCM, `lambda_outer.py`) → report, state continuum-vs-explicit limit.
- **λᵢ (inner-sphere, Nelsen 4-point, `reorg.py`)** — LOW CONFIDENCE: gas-phase
  on SMD geometries, single-conformer, negative values are a QC flag → DEMOTE to
  SI as qualitative, OR conformer-average before any main-text use. Do NOT make
  a kinetic/reversibility claim from it.
- **Disproportionation** — report as computed, conditional.

## Open requirements before submission (from PUBLICATION_READINESS.md)
- M2: re-confirm maintained code reproduces the frozen 44-couple benchmark.
- M7: one reproducible entry point per figure (config, seeds, env, raw-record
  location, content hashes); clean-env replay; path-relocation/restart test.
- Fig 6 depends on the running expand_catholyte campaign finishing + QC passing.
- Author metadata, data DOI, availability statement — user to supply.

## Reconcile with refactor plan
The pipeline-efficiency work (bag-of-tasks GPU queue, incremental tables) is the
engine behind "one reproducible entry point per result" (M7) and makes finishing
the Fig 6 catholyte expansion cheaper/faster. In a method-led paper the
efficiency gains are ALSO a legitimate methods contribution — report measured
GPU-utilization/throughput (but per M-acceptance: enumeration is not a
million-molecule QM benchmark; report real task costs, not a synthetic claim).
