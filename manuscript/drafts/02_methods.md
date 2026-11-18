> **Superseded (2026-11-18).** This per-section draft is kept for history only and is out of date:
> several numbers and claims in it were later corrected. The current text is
> `manuscript/manuscript_draft_v0.md` (journal version: `manuscript/jcim/`; SI:
> `manuscript/supporting_information.md`). Do not cite or copy from this file.

# 2. Computational Methods

We screen graftable redox motifs through a single reproducible driver
(`redox.workflows.motifs.commands`) that takes each registered motif through nine
stages: (1) build, (2) MLIP pre-optimization, (3) DFT geometry optimization with
implicit solvation per redox state, (4) refined single-point energies, (5)
thermochemistry, (6) derived redox properties, (7) reliability gating and
cross-point consistency, (8) multi-property ranking, and (9) an independent audit
of the exported tables. Each stage is addressed by content hash so that a change
in geometry, electronic protocol, or charge/multiplicity resolves a new record
rather than silently overwriting an old one. The level of theory and standard-state
conventions used throughout are collected in Table 1. Because this is a method-led
study, we state the protocol precisely here and defer the supporting derivations
(outer-sphere PCM operators, qRRHO projection, σ-aware domination proofs) to the
Supporting Information (SI).

## 2.1 Grafted-motif model

Each candidate is a redox motif attached to a single *4-methylbenzyl* cap that
serves as a Merrifield (poly(vinylbenzyl chloride)) backbone proxy. Grafting is
carried out as a balanced, local substitution onto 4-methylbenzyl chloride
(`Cc1ccc(CCl)cc1`); the ring methyl *para* to the benzylic linker stands in for the
polymer chain. Two reaction classes are recognized and charge/atom-balanced
explicitly (`structures/grafting.py`): N-quaternization of an azinium nitrogen
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

## 2.2 Conformer generation and MLIP pre-optimization

For each motif a charge-aware conformer ensemble is embedded with RDKit and
force-field ranked (`structures/candidates.py`). Each redox state is then pre-
optimized in the gas phase with the UMA machine-learned interatomic potential
(`backends/uma_current.py`), which is used purely as a fast, charge-/spin-aware
pre-optimizer to select the lowest-energy conformer and furnish a good DFT seed.
The backend is chosen from the *actual* reserved GPU architecture (deterministic
stock calculator on all current production cards); UMA is gas-phase only and never
supplies a solvated energy. By default a single best seed is carried into the DFT
stage, because the rigid aromatic redox cores collapse to one solvated minimum
(measured DFT seed spread < 5 meV in 13/14 reference states); multi-conformer
screening (`--nconf > 1`) and an optional ring–ring torsion scan are reserved for
genuinely floppy species and reduced dianions.

## 2.3 DFT geometry optimization and solvation

Every redox state is **optimized in solvent** on the SMD(acetonitrile) potential-
energy surface (geomeTRIC driver with grid-response gradients; `backends/dft.py`),
because a gas-phase MLIP seed cannot see solvation and charged/floppy species can
collapse without it. Each seed is first warm-started with a GFN2-xTB ALPB(MeCN)
pre-optimization, which also binds otherwise gas-unbound dianions. Geometries are
optimized at r²SCAN-D4/def2-SVP (anionic states def2-SVPD) in SMD. The SCF uses
RI-J density fitting (`conv_tol = 1 × 10⁻⁹`, `max_cycle = 200`) with an escalating-
stabilizer fallback (level shift/damping, then second-order Newton) for stubborn
open-shell cases; open-shell states are solved from several initial guesses
(minao/atom/Hückel) and the lowest converged solution is kept, with additional
orbital-removal occupation seeds for the Fe doublet of ferrocenium. A def2 ECP is
attached automatically for any element with Z ≥ 37. The three states of a motif are
independent and are optimized concurrently, one GPU each.

## 2.4 Refined electronic energies

Production energies are computed at a single uniform protocol,
**ωB97M-V/def2-TZVPD** with the VV10 nonlocal dispersion built into the functional,
SMD(acetonitrile), applied to *every* charge state and to the ferrocene reference
(`core/protocol.ACTIVE_SP`, `backends/sp.py`). Using one diffuse basis for all
charge states is deliberate: a ladder that scored neutrals with def2-TZVP and
anions with def2-TZVPD would place the diffuse freedom on one side of every
reduction and every disproportionation, so neither E° nor the solvation term would
get clean basis-set cancellation. Each energy is written to a protocol-addressed
record (`calcs/dft/<id>/<state>/sp/<hash>.json`) whose hash covers the coordinates,
charge, multiplicity, functional, basis, dispersion/NLC, solvent and solvent model,
density fitting, and convergence settings. A pre-existing optimizer energy is
*adopted* only if it proves it was produced at exactly this protocol (including
density-fitting provenance and the open-shell guess sweep); otherwise a fresh gas-
and SMD single point is computed. Both gas and solvated energies are retained, with
dG_solv = E_SMD − E_gas and a gas-HOMO/unbound-anion diagnostic.

## 2.5 Thermochemistry: a deliberate two-tier design

Free energies are G = E_SMD(active protocol) + G_thermal, and the two thermal terms
are maintained as two explicit tiers. The **pool-wide tier** is a GFN2-xTB RRHO
correction (298.15 K) evaluated on the DFT+SMD-optimized geometry
(`backends/dft._thermal_correction`): G_thermal = ZPE + H_thermal − TS. It is used
for ranking because RRHO corrections are nearly independent of the electronic-
structure method, so one cheap, uniform Hessian is applied identically to every
molecule and every state in the pool and cancels cleanly in the energy
*differences* we report. It is also the only safe option at scale here: the
installed gpu4pyscf analytic UKS Hessian inflates open-shell frequencies roughly
two-fold, corrupting ZPE for exactly the radical states that dominate this dataset.
The **lead-confirmation tier** is a native, derivative-checked solution DFT
thermochemistry (`thermochemistry/native.py`): an r²SCAN-D4/def2-SVP(D)/SMD Hessian
with full grid response, stationary-minimum and gradient acceptance checks, and a
spin-degeneracy term, with the nuclear partition evaluated by a certificate-locked
projected internal-mode qRRHO kernel (`qm/qrrho.py`). It is disk-/GPU-gated and
expensive, so it is run to *confirm* shortlisted leads rather than across the pool.
A state missing either an active energy or a valid thermal term makes its couple
INCOMPLETE; there is no electronic-only fallback, and a molecule and the ferrocene
reference must share one thermal definition or the couple is flagged
`thermal_definition_mismatch`.

## 2.6 Derived properties (kept as separate quantities)

Each electrochemical descriptor is computed and reported independently; they are
never combined.

- **Redox potential vs Fc.** For a one-electron couple O + e⁻ → R, dG = G(R) − G(O),
  E_abs = −dG, and E_vs_Fc = E_abs − E_ref, where the Fc/Fc⁺ reference is computed
  *live* at the active protocol with the same free-energy definition as the
  molecules, cancelling the systematic error of an absolute reference
  (`properties/potentials.py`; RT/F = 0.025693 V, Nernst slope 0.059160 V).
- **Usable capacity.** Specific capacity C_spec = nF/(3.6·MW) [mAh g⁻¹] on the
  accessible electron path, with electrons counted as *usable* only while their
  potentials lie inside the TBAPF₆/MeCN stability window ([−3.0, 2.0] V vs Fc,
  anolyte/catholyte divider 0 V); a second electron outside the window collapses to
  1e usable capacity.
- **Outer-sphere λ_o.** A linear-response continuum estimate by two methods
  (`properties/lambda_outer.py`): single-sphere Born/Marcus (SASA-equivalent
  radius) and a molecular-cavity nonequilibrium (Pekar) IEF-PCM using the real
  cavity and the charge-transfer surface potential, with ε_s = 37.5, ε_op = 1.806.
  Reported as the half-reaction (1-body) value; the continuum-vs-explicit gap is not
  closed.
- **Inner-sphere λ_i.** Nelsen four-point, gas-phase at the active protocol on the
  SMD geometries (`properties/reorg.py`). This is a current limitation: the gas
  energies on single-conformer SMD geometries make λ_i low-confidence (negative
  values and large symmetry-aware RMSD are QC flags), so λ_i is treated as
  qualitative/SI pending an all-SMD, conformer-matched surface. No kinetic or
  reversibility claim is drawn from it.
- **Disproportionation.** dG_disp = G(R^{q+1}) + G(R^{q−1}) − 2G(R^q) = F·(E_high −
  E_low), i.e. the reference-free two-wave spacing (`quality/stability.py`),
  reported as computed and conditional.

## 2.7 Reliability gating and quality control

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

## 2.8 Multi-property ranking

Ranking is performed on **σ-aware Pareto fronts with the anolyte and catholyte
pools kept separate** (`screening/pareto.py`, `screening/scorecard.py`). An
ambipolar molecule enters both pools with only its own side's electrons. The
objectives are the trustworthy axes only, oriented so that higher is better:
voltage (−E for anolytes, +E for catholytes), specific capacity, disproportionation
stability (only where an intermediate exists), and kinetics (−λ_i, QC-clean couples
only). Domination is σ-aware: A dominates B only if, over the objectives they share
(≥ 2), A is never worse beyond the combined per-axis noise k·√(σ_A² + σ_B²) and is
strictly better beyond noise on at least one axis, with each σ taken from the
candidate's own benchmark-derived uncertainty. A candidate missing any primary
objective is INCOMPLETE and cannot enter the front or dominate. Solubility and
synthetic-accessibility proxies, and the continuum λ_o/λ_het/λ_se, are carried as
annotations/tiebreakers, not objectives. There is deliberately **no scalar or
weighted-sum figure of merit**: such a score would require arbitrary weights, ignore
σ, and rescale as the candidate set changes; the front plus the raw axes are
reported instead.

## 2.9 Reproducibility and software

The pipeline is a command registry: `motifs.commands()` emits the exact argument
vector for every stage, and each run writes a `run_manifest.json` recording the
config and its hash, the selected-definition hash, host/PID/interpreter, backend,
GPU and CPU affinity, thread caps, xTB binary, conformer count, thermal definition,
and per-step return codes and timings; resume binds these inputs (append-only
registry growth is permitted, but a changed selected molecule fails). Numerical
reproducibility is anchored by five byte-exact portable replays (`examples/`, stdout
SHA-256 pinned) and a curated contract-test suite (744 passing). The qRRHO kernel is
certificate-locked (SHA-256-verified xTB thermodynamic library, isolated CPU
subprocess), as are the xTB kernel/runtime, the sampling snapshot, and the hash-
locked replay payloads; reformatting any of these invalidates the projected-thermal
certificate. All GPU work is submitted through a cluster-wide reservation gate
(reserve → run with `CUDA_VISIBLE_DEVICES` → release), with a per-motif file lock so
one writer owns each motif, `taskset` CPU pinning, OMP/MKL/OpenBLAS thread caps, and
a node-capacity/free-slot preflight; the gate is never bypassed. Final property
tables are accepted only after an independent arithmetic audit (zero discrepancies
required). An external experimental benchmark of the potential axis is maintained
(`quality/benchmark_panel.py`) and reported in Section [Results]; an end-to-end
re-run of the maintained code to re-confirm the frozen 39-endpoint benchmark is a
GPU/QM step deferred under the current no-QM constraint (reviewer item M2).

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
