# FINDINGS

Running log of what we've *learned* building this pipeline — the durable conclusions,
especially the negative results and the "do it this way, not that way" lessons. Task tracking
lives in `TODO.md`; this file is for insight.

Convention: each finding has the takeaway first, then the evidence.

---

## 0. The unifying principle (read this first)
**Energy *differences* with error cancellation compute reliably; *absolute* charged-species
solvation and phase transitions do not.** Everything below is a corollary.
- Reliable (cancellation): isodesmic/reference reactions, disproportionation (2R → R⁺+R⁻),
  reorganization energy λ (same molecule, two geometries), within-family *ranking*.
- Unreliable (absolute charged-species solvation / phase behavior): absolute redox potentials,
  the +2 pimer desolvation in dimerization, absolute solubility (needs the solid/lattice term),
  membrane crossover.

Design consequence: **build objectives on cancellation-friendly quantities; treat the rest as
labeled proxies or filters, never as trusted absolute numbers.**

---

## 1. Accuracy floor of implicit-solvent DFT redox potentials is ~0.5–0.9 V
Our Tier-1 pipeline: OROP experimental benchmark MAE **0.58 V** (n=36), vs OROP's *own* raw
implicit-DFT **0.43 V** on the same systems. This is the physics floor — the ~0.5 V "viologen
error" that started this was never anomalous; it is the normal accuracy of the method for
charged couples. **Ranking, not absolute potential, is the usable output.**

## 2. Ranking quality must be measured WITHIN charge class, not globally
Global Spearman (0.95) is inflated because charge classes separate (cations high, anions low),
so it "ranks" by charge, not chemistry. The honest metric is within-class:
- cations (+1): Spearman **0.885** — trustworthy for catholyte screening.
- anions (0/−1): Spearman 0.489 raw → see Finding 7.

## 3. gpu4pyscf's UKS (open-shell) analytic Hessian is BROKEN in this version
It inflates open-shell frequencies ~2× → corrupt ZPE (viologen radical cation: 14.7 eV vs
correct 6.4 eV). A contiguity monkeypatch made it *run* but the numbers were still garbage
(the bug is in the XC 2nd-derivative path, `_get_vxc_deriv2`). Closed-shell (RKS) is fine.
**Do not use GPU DFT Hessians for radicals.**

## 4. Thermal corrections: use GFN2-xTB RRHO (fast, correct, method-transferable)
Because RRHO thermal corrections are nearly method-independent AND the DFT Hessian is broken
for open shells (Finding 3), we compute G_thermal with xtb `--hess` (`dft._thermal_correction`).
Validated: viologen neu 5.13 eV vs ox1 5.21 eV (sane, consistent; matches the *correct* DFT
RKS ZPE). Effect on viologen: MAE 0.429 → **0.372**. Applied uniformly to the whole table.

## 5. NEGATIVE RESULT — released-counterion ion-pair scheme fails in implicit solvent
`2 MV²⁺·2PF6⁻` style neutral-assembly scheme (redox.ionpair). Even with a correct singlet
`mv_ip2`, waves/spacing are catastrophically wrong (spacing ~ +11 V vs exp −0.43). Cause:
forming a neutral ion pair from two well-solvated ions costs a large (~+6 eV), poorly-modeled
desolvation free energy that does NOT cancel between waves. **Do not resurrect continuum
ion-pairing.** Best physics viologen result = bare ions + thermal (MAE 0.372).

## 6. NEGATIVE RESULT — explicit microsolvation (4 MeCN) does not fix viologen
Cluster-continuum, all 3 states DFT+SMD-optimized. MAE 0.373 → 0.374 (no change). For a
*delocalized aromatic cation*, the error is not first-shell specific solvation, and MeCN is a
weak coordinating solvent. **Explicit solvation is not a general accuracy lever for redox
potentials (at least cations); it is expensive and was dropped.** (It may still matter for the
+2 dimer desolvation — Finding 9 — the one place it could help, but not pursued.)

## 7. KEY — the anion "ranking failure" was category errors, not a method failure
The anion Spearman collapse (0.489) was driven by ~5 molecules that **have no reversible
reduction**: CCl₄ and CH₂Br₂ (dissociative electron attachment — the radical anion breaks a
C–X bond; CCl₄ neutral→anion heavy-atom RMSD = 3.4 Å), an unbound radical anion (EA_gas < 0),
and a ring-opening anhydride. These were being scored against irreversible *peak* potentials —
a category error. Filtering to reversible couples restores anion ranking to **Spearman 0.95**
(MAE 0.53). A flow-battery anolyte *must* be reversible, so this filter is a required screening
criterion, not a convenience — and it doubles as a stability signal.

## 8. CAVEAT — a naive reversibility check has false positives (needs connectivity + solvated EA)
First implementation used whole-molecule heavy-atom RMSD < 0.8 Å + gas-phase EA > 0. It
FALSELY flagged:
- TEMPO (textbook *stable* radical) and functionalized anthraquinone as "dissociative" —
  because a floppy tail / ring pucker / methyl rotation inflates whole-molecule RMSD without
  any bond breaking.
- Anthraquinone's 2nd reduction (dianion) as "unbound" — because gas-phase dianions are almost
  always unbound (EA_gas < 0) even though they are perfectly bound and reversible in solution.
**Fix:** detect dissociation by **bond-connectivity change** between the two geometries (robust
to conformational flexibility), and test binding with **solvated EA (e_smd)**, not gas EA.
[status: fixing]

## 9. Stability = decomposition ΔG, but not all decomposition ΔGs are equally computable
- **Disproportionation** `2 R → R⁺ + R⁻`: atom- AND charge-symmetric → strong cancellation →
  RELIABLE. TEMPO• validated it (+221 kJ/mol, by far the most stable radical, as it must be).
  Exact from existing redox-state energies, no new calc. This is the trustworthy stability axis.
- **Dimerization** `2 R⁺ → dimer²⁺` (viologen pimer): needs the *absolute* solvation of a
  concentrated +2 → the continuum's weak spot → UNRELIABLE (±several tenths eV). Our value
  +0.49 eV (mildly stable) is consistent with MV⁺• being a persistent monomer, but the sign
  isn't firmly established. Report with a large error bar; lean on disproportionation.

## 10. Reorganization energy λ: standard, reliable, and validatable against D3TaLES
λ_i via the 4-point (Nelsen) scheme is IDENTICAL to D3TaLES's `ReorganizationCalc`
(verified in their source). Reliable via cancellation. Validation reference = D3TaLES computed
λ (28k molecules + SMILES): **hole column credible** (triphenylamine 0.096 eV ✓, phenothiazine
0.48 eV ✓ from its butterfly→planar flattening), **electron column noisy** (negatives, 10 eV
outliers — the same anion difficulty). Compute cross-points GAS-ONLY (`dft_smd(do_smd=False)`)
— inner-sphere λ needs only the gas energy, so skipping the SMD SCF ~halves the cost.

## 11. Campaign design: multi-objective Pareto, built correctly
Voltage + stability are necessary but not sufficient. Other axes and their computability:
- **capacity (n, MW): EXACT** (arithmetic, gated by the electrochemical window).
- **λ / kinetics: reliable** (cancellation).
- **solubility: relative only** (ΔG_solv proxy; absolute logS needs crystal lattice — hard).
- **SA / cost: cheap proxy** (heuristic, not a cost model).
- **crossover: descriptor only** (needs the membrane; not from single-molecule DFT).
Construct the front on the *trustworthy* axes; proxies are filters/annotations. Split into
anolyte/catholyte pools (redox potential is not monotonic-better). Use σ-aware domination
(don't rank within error bars). Pair the front with a physical figure of merit
(energy density ∝ n·V·solubility; $/kWh). Reversibility (Finding 7) is a hard pre-filter.

## 12. Explicit solvation is a Tier-2 *finisher*, never a screening tool
~10–20× the bare cost + shell sampling ⇒ hundreds of GPU-days for 1000 candidates. Use a
funnel: Tier-0 cheap ranking (thousands) → Tier-1 bare DFT+SMD+thermal (hundreds) → Tier-2
explicit solvation (dozens). "Apply broadly + physics-only + cheap": pick two.

## 13. Cross-project (metabolic ΔG pipeline) corroboration
The metabolic-reaction pipeline reaches MAE 13.5 kJ/mol ≈ **0.14 eV** with the same machinery —
6× better than our redox ~0.87 V — precisely because metabolic reactions are isodesmic/
group-conserving (cancellation), and it *maximizes* cancellation deliberately (cofactor-ring
swaps, truncation). Empirical proof that cancellation, not more physics, is the dominant lever,
and that same-charge/reference-reaction referencing is the right direction for redox.

## 14. System context: monomer as a proxy for a Merrifield-resin polymer
The real material is a **polymer** — a Merrifield resin (polystyrene backbone) with the
redox-active groups pendant via the chloromethyl/benzyl linker. We model the **monomer** because
the redox chemistry is *local* to each pendant group. Implications:
- The **benzyl / benzyloxy substituent IS the polymer tether**, not a real degree of freedom.
  Its floppiness in the isolated monomer (large inter-state RMSDs, conformer noise) is a
  **monomer-model artifact** — in the polymer the linker is constrained by the backbone. This
  vindicates using bond-connectivity (not RMSD) for the reversibility verdict (Finding 8) and
  argues for lightly restraining/ignoring the tether when it dominates conformer spread.
- **Solubility (Finding 11) does not transfer** monomer→polymer — polymer processability/swelling
  is a different property. Treat monomer "solubility" as low-confidence for this system.
- **Dimerization (Finding 9) stays relevant, possibly more so**: adjacent pendant radicals on
  the backbone are held at high local concentration, so inter-monomer π-dimerization is a real
  polymer failure mode the monomer model can only approximate.
- Redox potential, disproportionation stability, and λ are local and transfer reasonably.
- Strategy: **start with the monomer to establish we can model the local chemistry correctly**,
  then consider backbone/tether effects.

## 15. Stability axis validated vs experiment (disproportionation = wave spacing)
dG_disp = F*(E_high - E_low) is validated directly against experimental two-wave spacings
(redox.validate_stability): n=3 across families (MV 42, AQ 60, TEMPO 211 kJ/mol exp), MAE
**16 kJ/mol (0.17 eV)**, Spearman **1.00**, small +bias (we slightly over-stabilise). So the
disproportionation axis is trustworthy; its sigma is 0.17 eV. Broadening needs more molecules
with two measured MeCN waves (a compute task).

## 16. Selection engine assembled (scorecard -> Pareto shortlist)
- `scorecard.py`: unifies the 5 axis CSVs into one per-candidate row with sigma + trust per
  axis (config/scorecard_config.py). Gating: reversibility hard-filter; window gating of
  accessible n (WINDOW_V_VS_FC); anolyte/catholyte/ambipolar split; per-side potentials
  (a single mean is meaningless for ambipolar, e.g. TEMPO). Fc excluded (it's the reference).
- `pareto.py`: SIGMA-AWARE domination per pool (A dominates B only if never worse beyond
  combined noise, better on >=1) + transparent normalised figure of merit. Objectives =
  trustworthy axes only (voltage, capacity, stability, kinetics); SA/solubility are annotations.
- Lesson demonstrated: sigma-aware domination matters — phenothiazine has the top raw FoM
  (best voltage) but is DOMINATED by phenothiazine_parent (its voltage/lambda edges are within
  noise, parent's capacity is decisively higher). Raw scalar ranking would have misled.
- CAVEAT surfaced: TEMPO ranks top anolyte only via its APPROXIMATE, edge-of-window reduction
  couple (flagged) — a domain-judgement flag, not a trusted pick.

## 17. Operational notes
- xtb (GFN2) added to the env (setup_env.sh + check_env.py); provides all thermal corrections.
- The lambda cluster node is SHARED (another user's training jobs held GPUs 2/4 this session) —
  always check GPU ownership before launching; never assume "dedicated."
- ~20% of OROP DFT jobs hit SMD gradient non-convergence ("Nuclear gradients not converged") —
  needs an SCF-robustness patch before any large screening campaign.
