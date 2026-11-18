# Abstract, Conclusions, TOC graphic, and Associated Content

Target journal: Journal of Chemical Information and Modeling (JCIM).
Working title: *A validated, reproducible pipeline for screening graftable redox
motifs for non-aqueous organic redox flow batteries.*

All quantitative claims are drawn from SHARED_CONTEXT.md and manuscript_outline.md.
Items that could not be re-confirmed against the maintained code or that depend on
author-supplied metadata are marked **[TODO]**.

---

## Part 1 — Abstract (~190 words, single paragraph, JCIM style)

Selecting a redox-active motif to graft onto a preformed polymer backbone for a
non-aqueous organic redox flow battery (ORFB) requires redox properties computed
on the *exact* grafted species, gated for reliability, and anchored to experiment
— requirements that ad hoc density functional theory (DFT) screens rarely meet
together. We present an end-to-end, reproducible pipeline that builds exact
capped-graft models, computes redox potentials, usable capacity, and
reorganization components referenced to ferrocene (Fc) in acetonitrile, and
subjects every candidate to explicit quality-control gates before ranking. On an
unfitted retrospective benchmark the pipeline reproduces experimental redox
potentials with a mean absolute error of 0.247 V and a systematic +0.19 V bias
across 39 endpoints (30 molecules, 7 studies), so relative rankings are more
reliable than absolute values.
Applied to 91 candidates, the pipeline maps a strongly anolyte-dominated design
space (85 anolyte / 6 catholyte) and surfaces grafted 5,10-dihydrophenazine
two-electron catholytes that fill the sparse catholyte region, with DFT-derived
substituent design rules (electron-withdrawing N'-groups raise the potential,
ring methylation lowers it). Pairing these catholytes with pipeline anolytes
gives proposed full cells near 1.9 V. Code, byte-exact replays, and a contract-
test suite are released for independent reproduction.

---

## Part 2 — TOC graphic description (single figure)

A two-panel graphic. The left panel shows the grafting abstraction and pipeline
spine: a SMILES redox motif plus a 4-methylbenzyl (Merrifield-proxy) cap, taken
through a compact flow (build → UMA pre-optimization → ωB97M-V/def2-TZVPD with
SMD(acetonitrile) per redox state → quality gates → potential vs Fc). The right
panel shows a potential-versus-usable-capacity design map in which prior
candidates cluster on the anolyte side and newly surfaced grafted
5,10-dihydrophenazine two-electron catholytes (points colored by substituent)
land in the previously empty catholyte region, with an arrow indicating a
proposed ~1.9 V full cell. A single caption line states that a capped graft is a
model, not a synthesized material.

---

## Part 3 — Conclusions (~520 words)

We have presented a single, reproducible computational pipeline for screening
redox motifs intended for grafting onto preformed Merrifield-type polymers, and
demonstrated it on 91 candidates for non-aqueous organic redox flow batteries.
The contribution is the method rather than any one molecule: the pipeline builds
the *exact* grafted species as a 4-methylbenzyl-capped model, computes each redox
state at a consistent level of theory (ωB97M-V/def2-TZVPD with SMD in
acetonitrile), reports potentials against a pinned ferrocene reference together
with usable capacity and reorganization components, and admits a candidate to
ranking only after explicit quality-control gates on spin, geometry integrity,
state coverage, and cross-point consistency. Half-reactions and reorganization
terms are kept as separate reported quantities; candidates are ranked on
transparent, uncertainty-aware Pareto fronts computed separately for anolyte and
catholyte pools, with no single weighted-sum score.

The credibility of the applied results rests on an unfitted retrospective
benchmark: across 39 endpoints (30 molecules, 7 studies) the pipeline reproduces
experimental potentials with a mean absolute error of 0.247 V and a systematic
+0.19 V bias. Because the error
is largely a common offset, the pipeline ranks candidates more reliably than it
places absolute potentials. Against this anchor, the 91-candidate screen reveals
a strongly anolyte-dominated space (85 anolyte / 6 catholyte; 89 of 91 are
two-electron couples), and the pipeline's matched-hypothesis machinery surfaces
grafted 5,10-dihydrophenazine two-electron catholytes that populate the sparse
catholyte front. A matched single-substituent series yields interpretable design
rules — electron-withdrawing N'-substituents raise the potential versus Fc
(e.g., 2,2,2-trifluoroethyl to ca. +0.62 V) while ring methylation lowers it
(2,7-dimethyl to ca. +0.30 V). Pairing these catholytes with pipeline quinone and
imide anolytes gives proposed full cells near 1.9 V.

These results are bounded by explicit limitations. The proposed full cells are
proposals computed as the difference of two validated potentials; they are not
demonstrated capacity, cycling, crossover, or synthetic feasibility, and a
computed cell voltage does not establish any of these. A capped graft, a soluble
parent, and an oligomer are three distinct models, and the capped graft is not
the material. Inner-sphere reorganization (λ_i) is reported only qualitatively:
the present values rest on gas-phase energies evaluated on SMD geometries from
single conformers, so no kinetic or reversibility claim is drawn from them.
Occupied-polymer ensembles, neighboring-graft coupling, and synthetic routes are
out of scope for this study.

The outlook follows directly. A pipeline-efficiency refactor — a bag-of-tasks GPU
queue and incremental, content-hashed result tables — lowers the per-candidate
cost and makes larger prospective screens, and the completion of the catholyte
expansion, tractable. Lead candidates can then be confirmed on the pipeline's
second thermochemical tier, native DFT qRRHO, together with the in-progress
self-consistent, conformer-averaged treatment of nonequilibrium solvation and
reorganization. The pipeline, its byte-exact replays, and its contract-test suite
are released so that each reported number can be reproduced and extended
independently.

---

## Part 4 — Associated Content / Data and Code Availability

### Supporting Information

The Supporting Information (available free of charge at [TODO: DOI/URL]) contains:
the full method specification (level of theory, solvent model, standard state,
and ferrocene reference handling); the grafting abstraction and capped-graft
identity registry; the quality-control gate definitions and per-gate pass rates
(registered → state-complete → QC-clean → scored); the retrospective benchmark
panel — matched computed-versus-experimental potentials vs Fc in acetonitrile,
bootstrap intervals, and grouped-offset (endpoint-, molecule-, and study-held-out)
diagnostics; the complete 91-candidate scored table with per-candidate
uncertainties and family labels; the anolyte and catholyte Pareto fronts;
the matched 5,10-dihydrophenazine substituent series and derived design rules;
outer-sphere reorganization (λ_o) continuum estimates and the qualitative
inner-sphere (λ_i) treatment with its stated limitations; the proposed full-cell
pairings; and the numerical-replay, concurrency, and environment records.
JCIM requires a data-availability statement; see below.

### Code availability

Versioned code and curated metadata are maintained in the public project
repository (https://github.com/Ray16/RTFB-materials-research). The release provides the
`redox` command registry (unchanged public command names dispatching to the
documented package layout), five byte-exact portable replays under `examples/`
that reconstruct cached scientific results without rerunning DFT (verified via
`scripts/validation/check_refactor_oracles.py`), a scientific-contract test suite
of 773 passing tests (1 skipped; `pytest tests -q`, verified 2026-11-17), and a
certificate-locked projected-thermal qRRHO kernel whose byte-identity is required
for the replays.
Installation and a locked scientific environment are described by `setup_env.sh`
and `docs/INSTALLATION.md`. GPU work is gated through `gpu_reserve`; the released
tree runs reproduction and tests CPU-only and launches no DFT or Hessian jobs.

### Data availability

The raw quantum-chemical records (geometries, energies, frequencies, and state
metadata) for the candidate pool and the validation set reside in dedicated data
trees (~270 GB candidate campaign; ~52 GB benchmark/reference set) referenced by
the pipeline through a documented data root. A locally prepared reviewer bundle
reconstructs the primary benchmark, grafting, and finite-ensemble values offline
after relocation, with original cluster paths blocked; its scope, environment
records, and integrity checks are described in the Supporting Information. The
raw-record archive does not yet have a public archival DOI; final deposition and
author metadata remain to be supplied **[TODO: author to supply data DOI,
repository release tag, author list, contributions, funding, and competing-
interest statements]**.
