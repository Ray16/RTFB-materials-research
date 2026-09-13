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

## 10. Reorganization energy λ: standard, reliable, validated by REPRODUCING D3TaLES at their level
λ_i via the 4-point (Nelsen) scheme uses the same **formula** as D3TaLES's `ReorganizationCalc`,
but their **level of theory is different** and must be matched to compare numbers:
- **D3TaLES level = IP-tuned LC-ωHPBE / Def2SVP, gas** (Duke 2023; ω is tuned *per molecule* and
  stored in the dump's `omega` column). Their Def2SVP has **no diffuse functions** → their
  *electron* (anion) column is noisy (negatives, >3 eV outliers); their hole column is cleaner
  except for poor donors (e.g. quinone cations).
- **Our production level = ωB97M-V/def2-TZVP(D) on SMD-opt geoms, diffuse for anions** — a
  *different functional + basis + phase*, chosen because it's more appropriate for our reductive
  anolytes (diffuse functions are essential for anions).

**Do-it-correctly validation (`validate_reorg_worker_d3tales.py`):** reproducing their EXACT level
— `lc_wpbe` (= LC-ωHPBE) with *their* stored per-molecule ω + def2-svp, gas — matches their
reported λ closely on clean entries (naphthoquinone electron 0.5807 vs 0.5825; duroquinone hole
0.765 vs 0.777 & electron 0.546 vs 0.573). The disagreements are where **D3TaLES's own value is
the outlier** (naphthoquinone hole 0.168, 80TFSO electron 1.055 — implausible). So:
- The ~−0.24 eV offset our *production* λ shows vs D3TaLES is a **functional+basis difference, NOT
  a solvent effect** — it is identical in our gas and SMD calcs, and it collapses when we match
  their functional. (Superseded the earlier B3LYP/6-31G* "matched" run, which matched neither
  their level nor ours and has been removed.)
- For our screening we KEEP our level (diffuse-augmented, solvated) — it is the more correct
  treatment for anions; matching D3TaLES is only for the cross-check.

Compute cross-points GAS-ONLY (`dft_smd(do_smd=False)`) — inner-sphere λ needs only the gas
energy, so skipping the SMD SCF ~halves the cost.

**λ outliers = UNBOUND radical anions, not a bug (16/442, all thiosuccinimides).** A cluster of
flexible imides (maleimide–thiol adducts, motif `N–C(=O)–CH₂–CH(S–R)–C(=O)`) gave λ = 1.5–4.1 eV.
Diagnosis (QC now in the pipeline): the **vertical radical anion is UNBOUND** — the gas HOMO at the
neutral geometry is **positive** (+0.07 to +1.14 eV; ⟨S²⟩≈0.75, so not spin contamination), so the
`E_R_at_O` cross-point is a "neutral + free electron" energy, not a bound reduced state, and the
4-point λ is ill-defined. **No recompute recipe rescues them:** conformer-matching (frozen rotatable
dihedrals; RMSD 3.3→0.9 Å) leaves λ~4 (still uses the unbound point); SMD stabilizes <0.3 eV (not
enough to bind, λ stays ~3.6); even conformer-matched **and** SMD together stays 3.1–3.9 eV. **D3TaLES
gets the same** huge values independently (their gas def2-SVP), confirming it is a real low-EA
electronic-structure limit, not our error. Correct treatment = **flag + exclude**, not fabricate a λ
(these unbound anions are not viable reductions anyway). Scripts: `reorg_screen_unbound.py` (cheap
vertical-anion HOMO screen), `reorg_recompute_guarded.py` (conformer-matched + SMD attempt),
`finalize_reorg_qc.py` (writes `anion_unbound`/`reliable` into `comparison.csv`; parity figures plot
the 426 reliable, annotate the 16 excluded).

**Pipeline hardened so this can't pass silently again.** `dft_smd` now returns `gas_homo_eV`/
`anion_unbound` (HOMO>0) and `gas_s_squared`/`spin_contam`, and accepts a geomeTRIC `constraints`
file (frozen dihedrals → conformer-matched λ). `reorg.py` flags each couple `anion_unbound`,
`conformer_jump` (heavy-atom RMSD>0.4 Å), `negative_lambda`, or `negative_half`, recording
`rmsd_A`/`anion_homo_eV` in `results/reorganization.csv`.

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

## 18. ΔG_solv validated against EXPERIMENT — reliable for neutrals, fails only for concentrated charge
Direct experimental validation of our computed ΔG_solv = e_smd − e_gas (full SMD, CDS included),
comparing to two gold-standard databases with matched conventions (298 K, Ben-Naim 1 M→1 M,
kcal/mol): **FreeSolv** (water, open) and **MNSol** (acetonitrile, the set SMD was parameterized
against; used with MNSol's own M06-2X gas geometries — the textbook fixed-geometry SMD protocol).
Scripts: `validate_solvation.py`, `plot_solvation_validation.py` → `results/solvation_validation.csv`,
figures `solvation_parity.png` / `solvation_mae_summary.png`. n=100.
- **Neutrals (the screening-relevant property): MAE 0.79 kcal/mol** (water n=24 MAE 0.78 R² 0.93;
  MeCN n=7 MAE 0.84). Matches SMD's published ~1 kcal/mol accuracy → **our ΔG_solv is trustworthy
  for neutral solutes.** This answers the Sept-2 action item.
- **Ions decompose by CHARGE CHARACTER, not sign** — accuracy tracks charge delocalization:
  delocalized ions (protonated amines, carboxylates, phenolates) MAE ~3 kcal/mol with good ranking
  (anion R² 0.92); **small "hard" cations with concentrated charge** (protonated methanol, t-BuOH,
  H2S, acids) MAE **14.9** kcal/mol (systematically under-stabilized, up to +31). Continuum SMD
  cannot supply their strong first-shell H-bonding.
- **This is independent experimental confirmation of Finding 0**: continuum solvation is reliable
  for neutrals + delocalized/cancellation-friendly charge, and fails for concentrated charge — the
  same physics as the viologen dication (Finding 5). Solvent choice (water vs MeCN) does not change
  the neutral accuracy; charge concentration does.
- CAVEAT: experimental *neutral* solvation data in MeCN is genuinely scarce (MNSol has only 7);
  the water leg (n=24) carries the statistical weight for the neutral claim. MeCN's rich MNSol data
  is ionic (39 cations + 30 anions).

## 19. Starting-candidate potentials are the INTRINSIC solvated model — call out ion pairing
The starting-candidate batch (config/starting_candidates.py: ethyl viologen, PMDI, ammonium-NDI,
methoxy-quinones) is computed as **bare solvated molecular ions in implicit MeCN** (no explicit
PF6- / cation), i.e. the *intrinsic* molecular redox thermodynamics — consistent with Finding 5
(continuum ion-pairing refuted) and the right baseline. But that baseline is NOT the same as the
experimentally observed potential in a real PF6-/Li+ electrolyte, and the gap is largest for the
**second reduction of the imides**:
- **Ammonium-NDI:** literature (2024 JACS) shows cation (e.g. Li+) stabilization of NDI radical-
  anion/dianion **compresses the two reduction waves** (brings E2 up toward E1). Our `ndi_ammonium`
  keeps ONE tethered trimethylammonium on the non-graft imide N, so it *partially* captures that
  "cation stabilizes the anion" effect **intramolecularly** — a feature, but it means its computed
  E2 already includes some of the ion-pairing shift and should not be read as the ion-free value.
- **Reporting rule:** label these E° as "intrinsic solvated (no explicit counterion)"; treat the
  **first reduction / first oxidation as the more transferable number**, and the **second charging
  event (the concentrated -2 / +2) as ion-pairing-sensitive** with a larger error bar (Finding 0/18:
  concentrated charge is where continuum solvation is weakest). Ranking within the set still holds;
  absolute E2 vs a specific electrolyte needs explicit-ion or calibration work, not more continuum.
