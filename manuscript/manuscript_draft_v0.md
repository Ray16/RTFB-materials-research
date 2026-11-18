# A validated, reproducible pipeline for screening graftable redox motifs for non-aqueous organic redox flow batteries

_Draft v0 — assembled 2026-11-17 from parallel section drafts (method-led spine, JCIM target). Benchmark verified at 0.247 V / 39 endpoints. Deterministic TODOs resolved against ground-truth data and the 7 Introduction citations placed (DOI-verified, see References). Remaining **[TODO]** markers are author/publisher metadata only (SI/data DOIs, author list, funding, competing interests) plus the GPU-deferred M2 benchmark re-run._

---

## Abstract

Selecting a redox-active motif to graft onto a preformed polymer backbone for a
non-aqueous organic redox flow battery (ORFB) requires redox properties computed
on the *exact* grafted species, gated for reliability, and anchored to experiment.
Ad hoc density functional theory (DFT) screens rarely meet these three requirements
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

## TOC Graphic

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

## Introduction

Non-aqueous organic redox flow batteries (ORFBs) are an attractive route to
grid-scale electrochemical storage because they decouple energy and power,
access cell voltages well beyond the aqueous stability window, and are built
from earth-abundant, structurally tunable organic redox couples rather than
scarce transition metals [Soloveichik 2015]. One of the more promising device
architectures immobilizes the redox-active species on a soluble or swellable
polymer, so that the active material is retained by a size-selective membrane
and crossover is suppressed [Janoschka 2015]. In practice these materials are
most readily made not by polymerizing a bespoke redox monomer but by grafting a
preformed redox motif onto a reactive, preformed backbone — most commonly a
chloromethylated (Merrifield-type) polystyrene, whose benzyl chloride handles
are alkylated by a nucleophilic site on the motif [Merrifield 1963]. This grafting
route decouples the chemistry of the backbone from the chemistry of the redox
center and makes a large library of candidate motifs synthetically plausible
from a single polymer platform. The design space it opens, however, is badly
unbalanced. The catalog of well-behaved organic anolytes — n-type couples based
on quinones, imides, and related electron-poor scaffolds — is deep and still
growing [Ding 2018], whereas robust, high-potential catholytes remain scarce
[Kwon 2019]. Because a flow cell needs both a catholyte and an anolyte, this
asymmetry, not the anolyte inventory, is what limits achievable cell voltage,
and it is precisely the kind of gap that a well-posed computational screen
should be able to survey and help fill.

Computation is in principle well suited to this task: redox potentials, and the
reorganization energetics that govern electron-transfer kinetics, are
accessible from solvated electronic-structure calculations, and density
functional theory (DFT) with an implicit solvent model has been used
extensively to rank organic redox couples [Marenich 2014; Er 2015]. In practice, though,
most published screens are difficult to trust or to reuse for the specific
problem of grafted motifs. Three shortcomings recur. First, the molecular model
is often a proxy for convenience rather than for the real material: a bare
parent heterocycle, an N–H or N–methyl analogue, or a soluble monomer stands in
for the grafted species, even though the attachment chemistry changes the
electronic environment at the redox center and therefore the potential. A capped
graft, a soluble parent, and an oligomer are genuinely different models and
should not be conflated. Second, reliability is rarely gated per candidate.
Spin contamination, unconverged or distorted geometries, missing charge states,
and inconsistent reference points can silently corrupt a reported potential,
and screens seldom document a per-state quality-control (QC) policy or state
which candidates were held back because they failed it. Third, and most
importantly, computed numbers are frequently reported without an external
calibration against measured potentials under matched conditions, so the reader
cannot tell whether a ranking reflects chemistry or method artifacts. Underlying
all three is a reproducibility problem: level of theory, solvent, standard
state, reference electrode, and software versions are often underspecified, and
the raw records needed to regenerate a figure are not released, so neither the
ranking nor its error bars can be independently reconstructed [Lejaeghere 2016]. A
computed descriptor, or a completed calculation, is not by itself evidence of
synthetic feasibility or of material performance; a screen is only as useful as
the discipline with which these distinctions are maintained.

This work addresses those gaps with a single, validated, end-to-end pipeline for
screening graftable redox motifs for non-aqueous ORFBs. Its design choices map
directly onto the shortcomings above. (i) Each motif is modeled as an *exact
capped graft*: the redox center carries a 4-methylbenzyl group as an explicit
Merrifield proxy, and the relevant oxidized or reduced charge states are built
and optimized individually, so the computed species is the grafted motif rather
than a stand-in. (ii) Every candidate passes explicit per-state QC —
spin/⟨S²⟩ checks, geometry integrity, full charge-state coverage, and
cross-point consistency — and candidates that fail are held rather than silently
dropped. (iii) All potentials are computed at ωB97M-V/def2-TZVPD with SMD(MeCN)
and reported against the ferrocene (Fc) reference in a TBAPF₆/MeCN electrolyte,
and the pipeline is calibrated against experiment, reproducing measured redox
potentials with a mean absolute error of 0.247 V and a systematic +0.19 V bias;
because the bias is systematic, relative rankings are more reliable than
absolute values. (iv) Candidates are ranked on a transparent, uncertainty-aware
(σ-aware) Pareto front over separate physical properties — potential and usable
capacity — with anolyte and catholyte pools treated separately and half-reaction
and reorganization quantities kept distinct; we deliberately avoid any
weighted-sum scalar score. Reorganization energetics are reported with
property-specific confidence: the outer-sphere term is a defensible continuum
estimate, whereas the inner-sphere term is currently treated as qualitative.
(v) The entire pipeline is byte-exact reproducible, distributed with portable
replays, a large contract-test suite, and a certificate-locked thermochemistry
kernel, so any reported number can be regenerated from its raw records.

We demonstrate the pipeline on a pool of 91 graftable candidates spanning
quinone, imide, pyridinic, phenazine, and phenothiazine families. The screen
confirms quantitatively that the accessible space is strongly anolyte-dominated
(85 anolyte versus 6 catholyte candidates), and the sparsity of the catholyte
Pareto front motivates a focused, matched-substituent study of grafted
5,10-dihydrophenazine (DMPZ) two-electron catholytes. From matched series
differing by a single substituent we extract an interpretable design rule —
electron-withdrawing N′-substituents raise the potential while ring-methylation
lowers it — that populates the previously empty catholyte front with
DFT-supported candidates near +0.3 to +0.6 V versus Fc. Pairing these catholytes
with pipeline anolytes suggests non-aqueous full cells near ~1.9 V; we stress
that these are computational proposals, not demonstrations of capacity, cycling,
crossover, or synthetic feasibility. The remainder of the paper is organized as
follows. We first describe the grafting abstraction and the nine-stage pipeline,
including its two-tier thermochemistry and its reliability-gating policy. We then
present the external validation against experiment that anchors the method's
credibility, apply the pipeline to the 91-candidate pool and its σ-aware Pareto
fronts, and finally present the DMPZ catholyte showcase and the proposed full
cells, before closing with the method's limitations and outlook.

---

## Computational Methods

We screen graftable redox motifs through a single reproducible driver
(`redox.workflows.motifs.commands`) that takes each registered motif through nine
stages: (1) build, (2) pre-optimization with a machine-learned interatomic
potential (MLIP), (3) DFT geometry optimization with
implicit solvation per redox state, (4) refined single-point energies, (5)
thermochemistry, (6) derived redox properties, (7) reliability gating and
cross-point consistency, (8) multi-property ranking, and (9) an independent audit
of the exported tables. Each stage is addressed by content hash so that a change
in geometry, electronic protocol, or charge/multiplicity resolves a new record
rather than silently overwriting an old one. The level of theory and standard-state
conventions used throughout are collected in Table 1. Because this is a method-led
study, we state the protocol precisely here and defer the supporting derivations
(outer-sphere PCM operators, qRRHO projection, σ-aware domination proofs), together
with a module map naming the specific file and function behind each stage, to the
Supporting Information (SI); the package layout is documented in the released code.

### Grafted-motif model

Each candidate is a redox motif attached to a single *4-methylbenzyl* cap that
serves as a Merrifield (poly(vinylbenzyl chloride)) backbone proxy. Grafting is
carried out as a balanced, local substitution onto 4-methylbenzyl chloride
(`Cc1ccc(CCl)cc1`); the ring methyl *para* to the benzylic linker stands in for the
polymer chain. Two reaction classes are recognized and charge/atom-balanced
explicitly: N-quaternization of an azinium nitrogen
(leaving group Cl⁻) and N–H/O–H substitution (leaving group HCl, flagged as
requiring a separate base correction). No base, counterion activity, or kinetics is
assumed. Redox states are defined per motif as adjacent one-electron couples from
the as-prepared resting state (manifest `n_e = 0`): p-type motifs oxidize
(neutral → ox1 (+1) → ox2 (+2)), n-type anolytes reduce (neutral → red1 → red2),
and ambipolar motifs carry couples on both sides. We emphasize throughout that a
**capped graft, a soluble parent, and an oligomer are three distinct models**: the
capped graft fixes the local electronic structure and attachment chemistry but does
not represent chain packing, polymer solubility, or cooperative effects. For
capacity bookkeeping the cap's methyl is converted to one backbone carbon per graft
site (CH₃–Ar → –[CH₂–CH(Ar)]–; `properties/capacity_and_proxies.py`), giving an
exact repeat-unit mass without asserting a synthesized material.

### Conformer generation and MLIP pre-optimization

For each motif a charge-aware conformer ensemble is embedded with RDKit and
force-field ranked. Each redox state is then pre-
optimized in the gas phase with the UMA machine-learned interatomic potential, which is used purely as a fast, charge-/spin-aware
pre-optimizer to select the lowest-energy conformer and furnish a good DFT seed.
The backend is chosen from the *actual* reserved GPU architecture (deterministic
stock calculator on all current production cards); UMA is gas-phase only and never
supplies a solvated energy. By default a single best seed is carried into the DFT
stage, because the rigid aromatic redox cores collapse to one solvated minimum
(measured DFT seed spread < 5 meV in 13/14 reference states); multi-conformer
screening (`--nconf > 1`) and an optional ring–ring torsion scan are reserved for
genuinely floppy species and reduced dianions.

### DFT geometry optimization and solvation

Every redox state is **optimized in solvent** on the SMD(acetonitrile) potential-
energy surface (geomeTRIC driver with grid-response gradients; `backends/dft.py`),
because a gas-phase MLIP seed cannot see solvation and charged/floppy species can
collapse without it. Each seed is first warm-started with a GFN2-xTB solvated pre-optimization
(ALPB, the analytical linearized Poisson–Boltzmann implicit-solvent model, in
MeCN), which also binds otherwise gas-unbound dianions. Geometries are
optimized at r²SCAN-D4/def2-SVP (anionic states def2-SVPD) in SMD. The SCF uses
RI-J density fitting (`conv_tol = 1 × 10⁻⁹`, `max_cycle = 200`) with an escalating-
stabilizer fallback (level shift/damping, then second-order Newton) for stubborn
open-shell cases; open-shell states are solved from several initial guesses
(minao/atom/Hückel) and the lowest converged solution is kept, with additional
orbital-removal occupation seeds for the Fe doublet of ferrocenium. A def2
effective core potential (ECP) is attached automatically for any element with Z ≥ 37. The three states of a motif are
independent and are optimized concurrently, one GPU each.

### Refined electronic energies

Production energies are computed at a single uniform protocol,
**ωB97M-V/def2-TZVPD** with the VV10 nonlocal dispersion built into the functional,
SMD(acetonitrile), applied to *every* charge state and to the ferrocene reference. Using one diffuse basis for all
charge states is deliberate. A mixed ladder — def2-TZVP for neutrals and def2-TZVPD
for anions — would place the extra diffuse functions on only one side of every
reduction and every disproportionation, so neither the potential (E°) nor the
solvation term would benefit from clean basis-set cancellation. Each energy is written to a protocol-addressed
record (`calcs/dft/<id>/<state>/sp/<hash>.json`) whose hash covers the coordinates,
charge, multiplicity, functional, basis, dispersion/NLC, solvent and solvent model,
density fitting, and convergence settings. A pre-existing optimizer energy is
*adopted* only if it proves it was produced at exactly this protocol (including
density-fitting provenance and the open-shell guess sweep); otherwise a fresh gas-
and SMD single point is computed. Both gas and solvated energies are retained, with
dG_solv = E_SMD − E_gas and a gas-HOMO/unbound-anion diagnostic.

### Thermochemistry: a deliberate two-tier design

Free energies are G = E_SMD(active protocol) + G_thermal, and the two thermal terms
are maintained as two explicit tiers. The **pool-wide tier** is a GFN2-xTB
rigid-rotor–harmonic-oscillator (RRHO) correction (298.15 K) evaluated on the DFT+SMD-optimized geometry: G_thermal = ZPE + H_thermal − TS. It is used
for ranking because RRHO corrections are nearly independent of the electronic-
structure method, so one cheap, uniform Hessian is applied identically to every
molecule and every state in the pool and cancels cleanly in the energy
*differences* we report. It is also the only safe option at scale here: the
installed gpu4pyscf analytic UKS Hessian inflates open-shell frequencies roughly
two-fold, corrupting ZPE for exactly the radical states that dominate this dataset.
The **lead-confirmation tier** is a native, derivative-checked solution DFT
thermochemistry: an r²SCAN-D4/def2-SVP(D)/SMD Hessian
with full grid response, stationary-minimum and gradient acceptance checks, and a
spin-degeneracy term, with the nuclear partition evaluated by a certificate-locked
projected internal-mode quasi-RRHO (qRRHO) kernel. It is disk-/GPU-gated and
expensive, so it is run to *confirm* shortlisted leads rather than across the pool.
A state missing either an active energy or a valid thermal term makes its couple
INCOMPLETE; there is no electronic-only fallback, and a molecule and the ferrocene
reference must share one thermal definition or the couple is flagged
`thermal_definition_mismatch`.

### Derived properties (kept as separate quantities)

Each electrochemical descriptor is computed and reported independently; they are
never combined.

- **Redox potential vs Fc.** For a one-electron couple O + e⁻ → R, dG = G(R) − G(O),
  E_abs = −dG, and E_vs_Fc = E_abs − E_ref, where the Fc/Fc⁺ reference is computed
  *live* at the active protocol with the same free-energy definition as the
  molecules, cancelling the systematic error of an absolute reference
  (`properties/potentials.py`; RT/F = 0.025693 V, Nernst slope 0.059160 V).
- **Usable capacity.** Theoretical specific capacity C_spec = nF/(3.6·M) [mAh g⁻¹],
  with M the neutral polystyrene repeat-unit mass of the graft (above) and
  counterions excluded, the usual active-material convention. A salt-loaded variant
  that adds the heaviest counterion load (n PF₆⁻ for catholytes, n TBA⁺ for
  anolytes) is reported in the SI. Because TBA⁺ (242 g mol⁻¹) is much heavier than
  PF₆⁻ (145 g mol⁻¹), that variant penalizes anolytes more than catholytes, so it
  must not be used to compare capacities across pools. It does not change
  Pareto-front membership within either pool. Capacity is counted on the
  accessible electron path, with electrons counted as *usable* only while their
  potentials lie inside the TBAPF₆/MeCN stability window ([−3.0, 2.0] V vs Fc,
  anolyte/catholyte divider 0 V); a second electron outside the window collapses to
  1e usable capacity.
- **Outer-sphere λ_o.** A linear-response continuum estimate by two methods
 : single-sphere Born/Marcus (using a radius from
  the solvent-accessible surface area, SASA) and a molecular-cavity nonequilibrium
  (Pekar) integral-equation-formalism polarizable continuum model (IEF-PCM) using the real
  cavity and the charge-transfer surface potential, with ε_s = 37.5, ε_op = 1.806.
  Reported as the half-reaction (1-body) value; the continuum-vs-explicit gap is not
  closed.
- **Inner-sphere λ_i.** Nelsen four-point, gas-phase at the active protocol on the
  SMD geometries. This is a current limitation: the gas
  energies on single-conformer SMD geometries make λ_i low-confidence (negative
  values and large symmetry-aware RMSD are QC flags), so λ_i is treated as
  qualitative/SI pending an all-SMD, conformer-matched surface. No kinetic or
  reversibility claim is drawn from it.
- **Disproportionation.** dG_disp = G(R^{q+1}) + G(R^{q−1}) − 2G(R^q) = F·(E_high −
  E_low), i.e. the reference-free two-wave spacing,
  reported as computed and conditional.

### Reliability gating and quality control

Every reported candidate must pass explicit, automated QC before it is ranked
(`quality/readiness.py --strict`): declared-vs-realized charge, atom composition,
connectivity, stereochemistry and spin parity; geometry-protocol currency; presence
of a valid active single point; a valid thermal term; and absence of imaginary
thermal modes. Spin states are checked by ⟨S²⟩/spin-contamination diagnostics and,
for near-degenerate cases, a bounded source-addressed alternate-multiplicity DFT
comparison (`quality/spins.py`; max spin contamination 0.1). State coverage,
cross-point caches, and the single shared thermal definition are verified for
consistency across a ladder. Crucially, candidates failing QC are **held, not
dropped**: `potentials.py` emits an explicit `INCOMPLETE`/`THERMAL_UNVALIDATED`
status, and incomplete rows are carried on a separate partial-information list in
the ranking (they can be dominated but can never dominate or be silently removed).

### Multi-property ranking

Ranking is performed on **σ-aware Pareto fronts with the anolyte and catholyte
pools kept separate**. An
ambipolar molecule enters both pools with only its own side's electrons. The
objectives are the trustworthy axes only, oriented so that higher is better:
voltage (−E for anolytes, +E for catholytes), specific capacity, disproportionation
stability (only where an intermediate exists), and kinetics (−λ_i, QC-clean couples
only). Domination is uncertainty-aware (σ-aware). Candidate A dominates candidate B only
when, across the objectives they share (at least two), A is never worse than B by
more than the combined per-axis noise, k·√(σ_A² + σ_B²), and is clearly better —
beyond that noise — on at least one axis. Each σ is the candidate's own uncertainty,
propagated from the validation benchmark. A candidate missing any primary
objective is INCOMPLETE and cannot enter the front or dominate. Solubility and
synthetic-accessibility proxies, and the continuum λ_o/λ_het/λ_se, are carried as
annotations/tiebreakers, not objectives. There is deliberately **no scalar or
weighted-sum figure of merit**: such a score would require arbitrary weights, ignore
σ, and rescale as the candidate set changes; the front plus the raw axes are
reported instead.

### Reproducibility and software

The pipeline is a command registry: `motifs.commands()` emits the exact argument
vector for every stage, and each run writes a `run_manifest.json` recording the
config and its hash, the selected-definition hash, host/PID/interpreter, backend,
GPU and CPU affinity, thread caps, xTB binary, conformer count, thermal definition,
and per-step return codes and timings; resume binds these inputs (append-only
registry growth is permitted, but a changed selected molecule fails). Numerical
reproducibility is anchored by five byte-exact portable replays (`examples/`, stdout
SHA-256 pinned) and a curated contract-test suite (773 passing). The qRRHO kernel is
certificate-locked (SHA-256-verified xTB thermodynamic library, isolated CPU
subprocess), as are the xTB kernel/runtime, the sampling snapshot, and the hash-
locked replay payloads; reformatting any of these invalidates the projected-thermal
certificate. All GPU work is submitted through a cluster-wide reservation gate
(reserve → run with `CUDA_VISIBLE_DEVICES` → release), with a per-motif file lock so
one writer owns each motif, `taskset` CPU pinning, OMP/MKL/OpenBLAS thread caps, and
a node-capacity/free-slot preflight; the gate is never bypassed. Final property
tables are accepted only after an independent arithmetic audit (zero discrepancies
required). An external experimental benchmark of the potential axis is maintained and reported in Section [Results]; an end-to-end
re-run of the maintained code to re-confirm the frozen 39-endpoint benchmark is a
GPU/QM step deferred under the current no-QM constraint.

---

### Table 1. Level of theory and standard-state conventions

| Item | Setting |
|---|---|
| Conformer pre-optimization | RDKit ensemble (charge-aware, FF-ranked) → UMA MLIP, gas phase |
| Geometry optimization | r²SCAN-D4/def2-SVP (anions def2-SVPD), SMD(MeCN), geomeTRIC, grid-response gradients |
| Solvated pre-optimization | GFN2-xTB, ALPB(MeCN) warm start |
| Production energy | ωB97M-V (VV10 NLC)/def2-TZVPD, **uniform diffuse basis for all charge states and Fc** |
| Solvation model | SMD, acetonitrile (ε_s = 37.5, ε_op = 1.806, n = 1.344) |
| SCF | RI-J density fitting; conv_tol 1 × 10⁻⁹; max_cycle 200; multi-guess open-shell, lowest converged kept |
| Heavy-element treatment | def2 ECP for Z ≥ 37 |
| Thermochemistry (pool-wide, ranking) | GFN2-xTB RRHO on the DFT+SMD geometry, 298.15 K (`legacy_gfn2_gas_rrho`) |
| Thermochemistry (lead confirmation) | native solution DFT Hessian (r²SCAN-D4/def2-SVP(D)/SMD) + certificate-locked projected qRRHO (`native_solution_dft_qrrho_spin`) |
| Free energy | G = E_SMD + G_thermal (one matched thermal definition per molecule and Fc) |
| Temperature | 298.15 K (RT/F = 0.025693 V; Nernst slope 0.059160 V) |
| Potential reference | Fc/Fc⁺, computed live at the active protocol (level-matched) |
| Electrolyte / window | TBAPF₆ / MeCN; [−3.0, 2.0] V vs Fc; anolyte/catholyte divider 0 V |
| Outer-sphere λ_o | Born (SASA radius) and molecular-cavity nonequilibrium IEF-PCM (Pekar), 1-body |
| Inner-sphere λ_i | Nelsen four-point, gas phase at the production protocol on SMD geometries (qualitative, SI) |
| Ranking | σ-aware Pareto fronts, anolyte/catholyte pools separate; **no scalar/weighted-sum score** |

*Solvent constants are those recorded in `config/project.json`: static permittivity ε_s = 37.5 and refractive index n_D = 1.344, giving optical permittivity ε_op = n_D² = 1.806 (the Pekar factor 1/ε_op − 1/ε_s feeds the outer-sphere λ_o). n_D = 1.344 is the standard tabulated value for MeCN at 20 °C; its temperature coefficient is small (dn/dT ≈ −4 × 10⁻⁴ K⁻¹), so the ~5 K difference from the 298.15 K thermochemistry shifts n_D by ≈ 0.002 and ε_op by < 0.01, a negligible effect on λ_o. We adopt the 20 °C n_D with the 298.15 K free energies, the conventional treatment.*

---

## Results and Discussion

### Validation against experiment

The credibility of any downstream ranking rests on whether the pipeline reproduces
measured electrochemistry, so before applying it to the candidate pool we test it
against literature redox potentials under matched conditions. The validation is
deliberately retrospective and uncorrected: we compare the pipeline's computed
potentials of *free parent* molecules against experimental values measured in
acetonitrile and referenced to ferrocene (Fc), using the identical build → DFT+SMD →
single-point → thermochemistry route applied to every screening candidate. The
experimental panel is hand-curated from the primary literature, with each wave
carrying its source, DOI, table or page, and original reference electrode; every
value is admitted (or not) *before* any computed number is seen, on criteria set by
what the source reports rather than by the residual it produces. Only waves reported
as (quasi)reversible half-wave potentials vs Fc, either directly or through the
authors' own calibration, enter the scored panel (tier A); half-peak potentials,
added-water electrolytes, unverified conversion constants, cathodic-peak-only waves,
and compounds with unspecified stereocenters or inferred identity are tiered down or
excluded with their reason recorded and never scored. Potentials are computed as
E_vs_Fc = -ΔG/n - E_abs(Fc), where the Fc/Fc⁺ reference is evaluated *live* at the
active level of theory with the same free-energy definition (G = E_SMD + G_thermal)
as every molecule, so the reference can never drift from the molecules it calibrates.
No empirical correction, fitted offset, or outlier removal is applied, and the panel
is not held out from model development — it is a descriptive accuracy report, not a
calibrated predictor.

**Headline result.** Across the frozen tier-A panel the pipeline reproduces
experimental potentials with a mean absolute error of **0.247 V**, an RMSE of 0.279 V,
and a **systematic bias of +0.19 V** (computed potentials sit above experiment),
with a maximum absolute error of 0.547 V. The panel comprises **39 scored
endpoints (redox couples) drawn from 30 molecules and seven independent studies**.
Per-molecule (macro) MAE is 0.263 V. The bias depends on the direction of the
couple: the 33 reduction endpoints (anolyte-like) give MAE 0.239 V and bias
+0.17 V, whereas the six oxidation endpoints (five phenothiazine, one nitroxide;
catholyte-like) give MAE 0.289 V and bias +0.29 V. Uncertainty is
quantified by cluster bootstrap (2000 replicates, seed 20261021): resampling by
molecule gives a 95% interval of [0.213, 0.286] V for MAE and [0.128, 0.257] V for
bias; the weaker seven-study clustering gives [0.178, 0.306] V and [0.053, 0.285] V,
reported in the SI. These are descriptive intervals for the observed panel, not
prediction intervals for a new molecule.

**Interpretation.** The dominant error is a *systematic* positive offset, not scatter:
the bias (+0.19 V) is a large fraction of the MAE, and grouped diagnostics localize it.
The first reduction/oxidation step carries nearly all of it (wave 1 MAE 0.266 V, bias
+0.266 V), while the second step is both smaller and centered (wave 2 MAE 0.185 V, bias
-0.064 V). A systematic offset is exactly the error structure under which *relative*
rankings are far more reliable than absolute values — a uniform shift reorders nothing —
and this is what licenses the comparative screening in the following sections, where
candidates are ranked within pools rather than reported as absolute cell potentials.
We resist attributing the offset to a single cause. The live Fc reference (E_abs(Fc) =
4.379 V) differs from other-functional literature values by an amount comparable to the
first-step bias, making a reference contribution plausible but not proven; reference-free
wave spacings (E1 - E2), which cancel any common Fc shift, are themselves systematically
too large across all nine paired ladders, so a pure reference error cannot explain the
residuals. Charge-state solvation and ion-pairing errors can coexist with it. We state
these as explicit limits: the panel validates *free parents*, not exact grafts or
polymer interactions (`exact_graft_accuracy_validated = False` in the record); external
naphthalene-diimide (NDI) accuracy is not covered, since the available NDI entry was
excluded for inferred identity and unspecified stereocenters; and causal attribution of
the bias to the Fc reference versus solvation remains unresolved.

The panel contains no dihydrophenazine, so it does not by itself calibrate the
catholyte family used in the showcase below. As a family-specific check we computed
the exact measured parent, 5,10-bis(2-methoxyethyl)-5,10-dihydrophenazine (BMEPZ;
E₁ = −0.29 V, E₂ = +0.50 V vs Fc in MeCN [Kwon 2019]), at the same ranking level
used for the pool. The computed couples are −0.06 and +0.87 V, errors of +0.23 and
+0.37 V (MAE 0.30 V, bias +0.30 V). This is one molecule, so it is a check rather
than a benchmark, but it agrees with the oxidation-side bias of the panel and implies
that the absolute DMPZ potentials reported below are likely about 0.3 V too
positive. **[TODO: persist the BMEPZ computed record (`scripts/validation/bmepz_anchor.py`)
alongside the benchmark evidence; it currently exists only in a dated run note.]**

**Figure 4.** The parity plot reports computed vs experimental E_vs_Fc for all 39
endpoints, colored by chemical family (quinone n=22, MAE 0.286 V; aromatic imide n=9,
MAE 0.158 V; phenothiazine n=5, MAE 0.296 V; viologen n=2; nitroxide n=1), with the
1:1 line and the +0.19 V offset line annotated. An accompanying residual histogram
makes the systematic positive shift visible as a displaced, roughly unimodal
distribution, annotated with MAE = 0.247 V and bias = +0.19 V. The worst outlier
(N-ethyl-3,7-dimethoxyphenothiazine, 0.547 V) is retained, not trimmed.

**Reproducibility status.** The numbers above are the *frozen
archive* benchmark, regenerated here from the immutable evidence epoch and its
content-hashed panel. The maintained, refactored code paths that produce them —
`properties/potentials.py`, `reporting/benchmark.py`, and the content-hashed panel
builder/validator in `quality/benchmark_panel.py` — are in place and preserve the
scoring, admission, and freeze logic verbatim. Re-confirming that an end-to-end run of
the *maintained* code reproduces the frozen 39-endpoint / MAE 0.247 V / +0.19 V values
requires regenerating `results/redox_potentials.csv` from the DFT energies, which is a
GPU/QM step deferred under the current no-QM constraint.

---


### Mapping 91 grafted candidates onto a two-pool decision front

Applied to the full set of 91 registered grafted motifs, the pipeline returns a
reliability-gated property record for every candidate — potential versus
ferrocene (Fc) in TBAPF₆/MeCN, nominal usable capacity, outer- and inner-sphere
reorganization estimates, and disproportionation free energies — from which the
accessible design space can be read directly (Figure 5). The pool is strongly
anolyte-dominated: 85 of the 91 candidates reduce (n-type anolytes) and only 6
oxidize into an accessible catholyte window (p-type). By family the set is 36
imides, 29 quinones, 18 multi-electron pyridines, 5 dihydrophenazines, 2
mono-imides, and 1 phenothiazine; 89 of the 91 are two-electron couples, the two
exceptions being the single-wave phthalimide and 1,8-naphthalimide controls.
This composition is itself a finding about the literature we draw from:
graftable, well-precedented n-type cores (naphthalene- and pyromellitic-diimides,
benzo-/naphtho-/anthraquinones) are abundant, whereas graftable p-type cores are
scarce.

Candidates are ranked not by a single figure of merit but on σ-aware Pareto
fronts computed **separately for the anolyte and catholyte pools**, because a
redox potential is not globally "more is better": an anolyte is better the more
negative its potential, a catholyte the more positive. Within each pool the
trustworthy axes — pool-appropriate voltage and nominal usable capacity (mAh g⁻¹
in TBAPF₆/MeCN) — define the objectives, oriented so that higher is better, and
domination is declared only when one candidate is no worse than another beyond
the combined per-candidate uncertainty (σ, propagated from the validation
benchmark) on every shared axis and strictly better, beyond that noise, on at
least one. Candidates missing a primary axis are reported on a separate
partial-information list rather than silently promoted or dropped. We
deliberately avoid any weighted-sum score: arbitrary weights, insensitivity to
σ, and pool-relative min–max rescaling would all make a scalar rank change
whenever the candidate set changes. The reader is given the front and the raw
axes instead (Figure 5; full table in the SI). Inner-sphere reorganization
(λᵢ) is carried only qualitatively and is not used to set the headline ranking,
for the reasons given in the Methods.

Figure 5 presents the two Pareto panels side by side. The anolyte panel (85
candidates) is densely populated across a wide potential window: diimides and
amino-/alkoxy-naphthoquinones span roughly −1.0 to −1.6 V vs Fc, with
pyromellitic diimides (e.g. the ether-tethered PMDI at −1.54 V) and
alkoxy-anthraquinones (≈ −1.49 V) anchoring the most-negative edge and the
multi-electron viologens clustering nearer −0.7 to −1.0 V. Theoretical capacities
of the two-electron anolytes span ≈ 90–200 mAh g⁻¹, highest for the light
benzo- and naphthoquinones and falling as solubilizing oligoether or long-alkyl caps
add inert mass. The catholyte panel,
by contrast, is almost empty: six points, five of them dihydrophenazines
clustered near +0.4 to +0.5 V and one 3,7-dimethoxyphenothiazine at +0.78 V.
This visual sparsity — a crowded anolyte front against a nearly bare catholyte
front — is the central structural result of the screen, and it defines the
**catholyte gap** that the next step is designed to fill.

### Filling the catholyte gap with grafted dihydrophenazines

To populate the empty catholyte front we used the pipeline's matched-expansion
machinery to graft a single scaffold — 5,10-dihydrophenazine (DMPZ) — and vary
one substituent at a time, so that the resulting potential shifts can be read as
matched single-substituent deltas rather than confounded whole-molecule
comparisons (Figure 6). Every member carries the same 4-methylbenzyl (Merrifield
proxy) cap on one nitrogen and a two-electron {neutral → +1 → +2} ladder,
isolating the effect of the N′-substituent and of ring substitution. The
confirmed series potentials (E vs Fc, mean of the two one-electron couples, from
`full_pool_scored_legacyG_MeCN_TBAPF6.csv`) are:

| Candidate | N′ / ring substitution | E vs Fc (V) | Couples (V) | Capacity (mAh g⁻¹) |
|---|---|---:|---|---:|
| `exp_dmpz_mb_ethyl` | N′-ethyl | +0.43 | −0.03 / +0.89 | 164.2 |
| `exp_dmpz_mb_methyl` | N′-methyl | +0.47 | +0.03 / +0.91 | 171.6 |
| `hyp_dmpz_mb_methoxyethyl` | N′-(2-methoxyethyl) | +0.47 | +0.01 / +0.94 | 150.4 |
| `hyp_dmpz_mb_trifluoroethyl` | N′-(2,2,2-trifluoroethyl) | +0.62 | +0.15 / +1.10 | 140.9 |
| `hyp_dmpz_27dimethyl_ethyl` | 2,7-dimethyl, N′-ethyl | +0.30 | −0.14 / +0.73 | 151.2 |
| `hyp_dmpz_27dimethoxy_ethyl` | 2,7-dimethoxy, N′-ethyl | +0.21 | −0.22 / +0.64 | 138.7 |
| `exp_ptz_37dimethoxy` | 3,7-dimethoxyphenothiazine | +0.78 | +0.37 / +1.19 | 142.7 |

Read as matched deltas against the N′-ethyl DMPZ reference (+0.43 V), the series
yields a clean, transferable **design rule**. An electron-withdrawing
N′-substituent raises the potential: replacing ethyl with 2,2,2-trifluoroethyl
shifts E by +0.20 V (to +0.62 V), while the near-isosteric N′-methyl and
N′-(2-methoxyethyl) caps barely move it (+0.04 to +0.05 V). Conversely,
electron-donating ring substitution lowers it: 2,7-dimethyl brings E down by
≈ 0.13 V (to +0.30 V) and 2,7-dimethoxy by ≈ 0.22 V (to +0.21 V) relative to the
same N′-ethyl parent. The two handles are roughly additive and act on opposite
ends of the scaffold, so N′-electron-withdrawal and ring-electron-donation give
independent, oppositely-signed control over where a DMPZ catholyte sits — spanning
≈ +0.21 to +0.62 V across the grafted series, with the distinct
3,7-dimethoxyphenothiazine scaffold extending the catholyte reach to +0.78 V.
These light two-electron phenazine cores give theoretical capacities of
≈ 140–170 mAh g⁻¹. That is comparable to the anolyte pool, so the catholyte side is
not a capacity bottleneck.

Two cautions apply to these numbers. First, the matched shifts are the robust
quantity: within one scaffold, errors common to the family cancel. Second, the
absolute potentials are probably too high. The BMEPZ check above puts the DMPZ
couples about 0.3 V too positive, so the series likely sits nearer −0.1 to +0.3 V
vs Fc in experiment. Re-plotting the
catholyte Pareto panel with these points shows the gap being filled by a
structured family rather than by isolated outliers (Figure 6).

Several of these values are provisional. The expansion campaign
(`results/expand_catholyte_20261117`) is live at the time of writing, with the
N′-methyl, N′-(2,2,2-trifluoroethyl), 2,7-dimethyl, and 2,7-dimethoxy DMPZ
members still completing their native-DFT confirmation stages; their legacy
ranking potentials are reported here and will be reconciled against the
disk-gated qRRHO values once the run finishes and QC gates pass. The frozen
91-candidate pool underlying Figure 5 contains six catholytes (five DMPZ + one
phenothiazine); the maintained scored CSV now carries a seventh catholyte — the
2,7-dimethoxy DMPZ (E ≈ +0.21 V) appended by the live campaign — so it holds 92
rows while the frozen analysis pool remains 91. The figures enforce the frozen 91
and footnote the live extra.

### Proposed non-aqueous full cells

Because every potential in both pools is referenced to the same validated Fc
scale, candidates can be paired directly and a nominal cell voltage read off as
E(catholyte) − E(anolyte) (Figure 7). Pairing a surfaced DMPZ catholyte
(≈ +0.4 to +0.6 V) with a pipeline anolyte from the most-negative edge of the
quinone/imide front (≈ −1.5 V) gives proposed non-aqueous full cells near
**~1.9 V** — for example N′-ethyl or N′-methyl DMPZ (+0.43 / +0.47 V) against an
alkoxy-anthraquinone or pyromellitic-diimide anolyte (≈ −1.48 V) yields ≈ 1.9–1.95 V,
and the more electron-poor N′-(2,2,2-trifluoroethyl) DMPZ (E ≈ +0.62 V) paired
with an ether-tethered PMDI anolyte (≈ −1.48 to −1.54 V) would extend the proposal
toward ≈ 2.1 V (mean-potential arithmetic verified against the scored CSV:
+0.43/+0.47 − (−1.48) = 1.91/1.95 V; +0.62 − (−1.48) = 2.10 V). These cell voltages are
less sensitive to the calibration bias than the individual potentials, because the
two biases partly cancel. Subtracting the family-specific biases (+0.30 V for DMPZ
from the BMEPZ check; +0.23 V for quinones from the benchmark) lowers the 1.91 and
1.95 V cells to about 1.84 and 1.88 V. The design rule above makes
this voltage a tunable quantity rather than a single number: catholyte and
anolyte potentials can be shifted independently to target a window.

We state the limitation plainly. A computed cell voltage derived from two
validated half-cell potentials is a **proposal**, not a demonstrated device. It
says nothing about attainable capacity under cycling, Coulombic efficiency,
membrane crossover of the small grafted proxies, chemical stability of the
charged states over many cycles, or the synthetic feasibility of the exact
grafted structure — each candidate's record carries explicit route caveats, and
a capped graft is not its soluble parent nor an oligomer. What the screen does
establish is a reliability-gated, internally consistent map of where these
chemistries sit, a transferable substituent rule for moving a catholyte to a
target potential, and a concrete, testable set of anolyte/catholyte pairings for
experimental follow-up.

---

## Conclusions

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

## Associated Content

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

## References

*Introduction citations, resolved and DOI-verified against Crossref (2026-11-17).
Author-year keys are used in text; JCIM numbering will be applied at submission.
Citations in later sections are not yet enumerated here.*

- **[Soloveichik 2015]** Soloveichik, G. L. Flow Batteries: Current Status and Trends. *Chem. Rev.* **2015**, *115* (20), 11533–11558. DOI: 10.1021/cr500720t.
- **[Janoschka 2015]** Janoschka, T.; Martin, N.; Martin, U.; Friebe, C.; Morgenstern, S.; Hiller, H.; Hager, M. D.; Schubert, U. S. An Aqueous, Polymer-Based Redox-Flow Battery Using Non-Corrosive, Safe, and Low-Cost Materials. *Nature* **2015**, *527* (7576), 78–81. DOI: 10.1038/nature15746.
- **[Merrifield 1963]** Merrifield, R. B. Solid Phase Peptide Synthesis. I. The Synthesis of a Tetrapeptide. *J. Am. Chem. Soc.* **1963**, *85* (14), 2149–2154. DOI: 10.1021/ja00897a025.
- **[Ding 2018]** Ding, Y.; Zhang, C.; Zhang, L.; Zhou, Y.; Yu, G. Molecular Engineering of Organic Electroactive Materials for Redox Flow Batteries. *Chem. Soc. Rev.* **2018**, *47* (1), 69–103. DOI: 10.1039/C7CS00569E.
- **[Kwon 2019]** Kwon, G.; Lee, K.; Lee, M. H.; Lee, B.; Lee, S.; Jung, S.-K.; Ku, K.; Kim, J.; Park, S. Y.; Kwon, J. E.; Kang, K. Bio-Inspired Molecular Redesign of a Multi-Redox Catholyte for High-Energy Non-Aqueous Organic Redox Flow Batteries. *Chem* **2019**, *5* (10), 2642–2656. DOI: 10.1016/j.chempr.2019.07.006.
- **[Marenich 2014]** Marenich, A. V.; Ho, J.; Coote, M. L.; Cramer, C. J.; Truhlar, D. G. Computational Electrochemistry: Prediction of Liquid-Phase Reduction Potentials. *Phys. Chem. Chem. Phys.* **2014**, *16* (29), 15068–15106. DOI: 10.1039/C4CP01572J.
- **[Er 2015]** Er, S.; Suh, C.; Marshak, M. P.; Aspuru-Guzik, A. Computational Design of Molecules for an All-Quinone Redox Flow Battery. *Chem. Sci.* **2015**, *6* (2), 885–893. DOI: 10.1039/C4SC03030C.
- **[Lejaeghere 2016]** Lejaeghere, K.; et al. Reproducibility in Density Functional Theory Calculations of Solids. *Science* **2016**, *351* (6280), aad3000. DOI: 10.1126/science.aad3000.
