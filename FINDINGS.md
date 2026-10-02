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
**Historical result; superseded by Findings 22–23 for the production protocol.**
Our Tier-1 pipeline: OROP experimental benchmark MAE **0.58 V** (n=36), vs OROP's *own* raw
implicit-DFT **0.43 V** on the same systems. This is the physics floor — the ~0.5 V "viologen
error" that started this was never anomalous; it is the normal accuracy of the method for
charged couples. **Ranking, not absolute potential, is the usable output.**

## 2. Ranking quality must be measured WITHIN charge class, not globally
**Principle retained; numerical values superseded by Findings 22–23.**
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
`2 MV²⁺·2PF6⁻` style neutral-assembly scheme (`archive/src/redox/ionpair.py`). Even with a correct singlet
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
**Historical subset; the integrity lesson is retained, but current full-set metrics are in
Findings 22–23.**
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
- **Our production level = ωB97M-V/def2-TZVPD on SMD-opt geoms (uniform diffuse basis since #22; previously diffuse for anions only)** — a
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
**SUPERSEDED by #24:** the 3 spacings below were uncited (MV from an aqueous couple, TEMPO's
reduction approximate). Against 10 sourced tier-A spacings ΔG_disp is OVER-estimated by
~24 kJ/mol (all errors positive); σ_disp is now 0.29 eV.

dG_disp = F*(E_high - E_low) is validated directly against experimental two-wave spacings
(redox.validation.stability): n=3 across families (MV 42, AQ 60, TEMPO 211 kJ/mol exp), MAE
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

## 20. QC-AUDIT — lambda quality flags were computed but never reached the scorecard
An external review of the workflow raised eight issues; all eight were reproduced against the
code. The three load-bearing ones and what was done:

**(a) Flags computed, then dropped.** `redox.properties.reorg.lambda_for_couple` emits `flag`, `rmsd_A`
and `anion_homo_eV`, but `redox.screening.scorecard` built its lambda map from `lambda_i_eV` alone and
never read `flag`. The shipped `results/reorganization.csv` also predated those columns
(`rmsd_A`/`anion_homo_eV` absent). Re-running the current QC over the candidate set flags
**14 of 24 couples** — both ethylviologen couples, both ndi_ammonium, both mophquinone, and
`pmdi/red1->red2`. So the previously reported lambda means and the Pareto front were resting
on contaminated values. FIXED: flags now propagate; a flagged couple is excluded from the
mean, and the scorecard carries `lambda_qc` / `n_lambda_flagged` / `lambda_flags` plus an
unfiltered column for reference. Three of six starting candidates currently have NO QC-clean
lambda (`all_flagged`) and are therefore withheld from the Pareto plane rather than ranked.

**(b) The cache silently pinned the old protocol — the anion screen was DEAD CODE.**
`_cross_energy_gas` accepted any cache containing `e_gas_eV`. All 48 candidate cross-point
caches were written before the gas-HOMO/unbound-anion diagnostics existed, so
`bool(d.get("anion_unbound"))` evaluated `bool(None)` -> `False` for **every** couple: the
unbound-anion screen has never once fired on this set, and the cache would never regenerate
because the `e_gas_eV` test kept passing. The 14 flags above are therefore a LOWER BOUND.
FIXED: cache entries now carry `cache_schema`, a geometry content hash, charge/mult and the
required diagnostic fields; `_cache_ok` rejects anything else and `stale_cross_points()` sizes
a regeneration before running it. This matters most for the quinone dianions.

**(c) The large RMSDs are TORSIONAL, not core distortion — the data is sound.** Cartesian
heavy-atom RMSD reaches 2.6 A, which looks alarming, but in internal coordinates every flagged
couple shows **max bond-length change <= 0.042 A** (textbook inner-sphere distortion) against
**torsion changes of 70-180 deg**. So the DFT geometries/energies are fine; what is wrong is
the lambda DEFINITION (independent global minimization charges a tether/ring rotation to
lambda) and a whole-molecule RMSD test too blunt to tell the two apart. No electronic
structure needs recomputing. Several torsion changes are exactly 180 deg, which for a 2-fold
symmetric aryl ring is a symmetry-equivalent flip (physically the same conformer) — meaning
the flag partly OVER-triggers. Outstanding: automorphism-aware RMSD, then core-only vs
tether-only atom maps and a conformer-matched core lambda for ranking.

**Two traps found while verifying, not in the original review:**
- `build_candidates.py` embeds the pre-canonicalization RDKit mol but writes the CANONICAL
  SMILES to the manifest, so **manifest-SMILES atom indices do not map to the stored xyz atom
  order**. Any substructure match used to index into a geometry is silently scrambled. (Atom
  ordering IS consistent BETWEEN states of a molecule — element sequence and bond graph match
  exactly — so cross-state RMSD itself is valid.)
- `finalize_after_dft.sh` and `run_uma.sh` `source ~/miniforge3/...`; on lambda5 `$HOME` is
  `/homes/rzhu`, a DIFFERENT filesystem from the NFS repo home, so they cannot find the
  `redox` env at all. Use `/nfs/lambda_stor_01/homes/rzhu/miniforge3`.

**Lower-severity items confirmed:** spin multiplicity is taken from UMA with no DFT-level
confirmation (real weakness, but bit nothing here — all six candidates have spin gaps
0.41-1.32 eV, far above the 0.217 eV degeneracy threshold, so 0 of 6 are ambiguous);
`lambda_eV` was named as if total but is inner-sphere only (renamed `lambda_i_eV`); the
"before grafting" column is a methyl-capped analogue, not the literal precursor (relabelled
"minimal capped analogue"); Ertl SA does not measure graftability — **all 6 of 6** grafted
structures score LOWER than their capped analogue despite ~90 more Da, so dSA-on-grafting is
meaningless as a synthetic-cost signal and a route-level graftability score is still needed;
and `finalize_after_dft.sh` omits reorg/reversibility/capacity/scorecard/Pareto, which is
exactly how stale tables coexisted with newer QC code.

## 21. Batch-2 scope: dtbc_phenol is reported but NOT ranked (same rule as the metal-oxo rows)
`config/merrifield_multielectron.py` adds five grafted candidates from
`merrifield_multielectron_smiles.xlsx`. Two required repair and one is mechanistically out of
scope for a Marcus/Nelsen treatment:
- **dtbc_phenol** (`rankable=False`): the sheet grafts one catechol -OH as a benzyl ether then
  claims the 2e catechol/o-quinone couple — impossible, that couple needs BOTH oxygens. What
  remains is a 1e phenol oxidation, and even that is an **EC process** (ArOH-+ has pKa <~ -2 in
  MeCN, so electron transfer is followed by O-H deprotonation). The Nelsen 4-point lambda_i
  assumes one molecule on two adiabatic surfaces with NO bond made or broken; O-H cleavage
  violates that exactly as M=O formation does for the excluded metal-oxo rows. Computed and
  reported for completeness; excluded from the Pareto kinetics axis for consistency.
- **aq_benzylamino**: the sheet's C9-NH2 tether is destroyed by the anthracene->anthraquinone
  oxidation (C9 becomes the carbonyl), so it is modelled at C2. Checked for the obvious PCET
  risk and it is CLEAN — N-H stays 1.004/1.010/1.013 A across neu/red1/red2 with the nearest
  acceptor 4.7-6.0 A away, i.e. no intramolecular proton transfer. Residual risk is
  intermolecular (a neighbouring AQ(2-) deprotonating the N-H), which is bimolecular and
  invisible to a single-molecule model.
- The five transition-metal rows are excluded; see the module docstring. Their "resin-attached"
  SMILES are placeholders (Fe- and Mn-THPP both give the H2THPP FREE BASE, W gives WO2(acac)2,
  Re gives methyltrioxorhenium), SMD has no element-specific parameters for them (the MNSOL CDS
  term returns the same value to within 0.002 eV for Fe/Mn/Co/Mo/W on a fixed test geometry),
  and the Mn/Mo/W/Re couples are catalytic oxo transfer, not outer-sphere.

**Latent bug fixed in passing:** `redox.qm.dft` built every `gto.M` without an `ecp=`, so any
element from Rb (Z>=37) up would have been run ALL-ELECTRON in a valence-only def2 basis
(W: 74 electrons in the 40 AOs meant for 14) — it converges and returns a number rather than
erroring. `_ecp_for()` now attaches the matching def2 ECP only when a heavy element is present,
so all-light systems (every current molecule) are bit-for-bit unchanged.

## 22. QC-AUDIT II — lessons converted into enforced data contracts (2026-09-29)
A second external review plus our own audit of the computable axes. Every item below is now
enforced in code and covered by `tests/test_contracts.py` (21 tests, step T of the finalize
chain), not just documented.

**Silent-failure bugs that were changing results**
- **Missing data was scored as a negative.** The integrity gate looked up a verdict per couple
  and treated "no row" as "not reversible": the four batch-2 quinone/viologen candidates
  (textbook reversible 2e- anolytes) were REJECTED only because `reversibility.csv` predated
  them. Missing now = INCOMPLETE everywhere (scorecard, Pareto, E°, dG_disp, integrity).
- **Fc/Fc+ reference used a different free-energy definition.** Molecules used
  G = E_smd + g_thermal; the Fc reference used E_smd only -> every E° shifted by -79 mV
  (ferrocene read -0.079 V vs itself). The setter script regex-patched `electrolyte.py`,
  where the value no longer lived, and printed "[patched]" while changing nothing. The
  reference is now computed LIVE from the ferrocene states at the active protocol.
- **Mixed basis in every redox ladder.** Neutral/cation states were scored at def2-TZVP and
  anions at def2-TZVPD, so the diffuse freedom sat on one side of each reduction and each
  disproportionation. Uniform def2-TZVPD (`redox.core.protocol.ACTIVE_SP`, incl. Fc) shifts every
  first reduction (0 -> -1) by **-0.09 V (range -0.05 to -0.14 V, n=17)**, oxidations
  (+1 -> 0) by +0.016 V, second reductions 0 (both sides were already diffuse), Fc by
  +0.2 mV. So n-type E1 AND dG_disp (via E1 - E2) carried a systematic ~0.09 V error.
- **The OROP benchmark was not measuring the production protocol.** 134/170 system pairs had
  no thermal term (read as 0), and one pair had it on one side only (an error of eV). Backfilled
  (`scripts/validation/orop/orop_backfill_thermal.py`), pairs without thermal are now excluded,
  and the benchmark is re-scored at the same uniform protocol (`redox.qm.sp --orop`).
- **PCM lambda_o had an inverted sign hidden by abs().** In PySCF q = K^-1 R v with
  R = -f(eps), so the Pekar term is 1/2 dV.(q_op - q_s) >= 0; the code computed the negative
  and took abs(). Magnitudes were right, but any real operator error would have been invisible.
  Fixed; a negative value now raises; a genuine limiting-case test (unit point charge in a
  single-atom cavity) reproduces analytic Born to <1e-5. The previous "sphere_sanity_check"
  never ran PCM.
- **UMA spin gaps are not reliable enough to set confidence.** aq_benzyloxy dianion S-T gap:
  UMA 0.19 eV (flagged near-degenerate) vs DFT 0.81 eV (singlet confirmed). Spin confidence
  now uses a DFT check (`calcs/spincheck/`) when present.

**Definitions that were scientifically loose**
- lambda_total added the 4-point sum (lambda_O + lambda_R, the inner term of a SELF-EXCHANGE
  pair) to a 1-body (electrochemical) lambda_o — neither convention. Now reported explicitly:
  lambda_O, lambda_R, their sum, lambda_het = (lambda_O+lambda_R)/2 + lambda_o,1 (electrode ET)
  and lambda_se(d) = lambda_O+lambda_R + 2 lambda_o,1 (1 - a/d) (two-sphere self-exchange).
- "reversible" over-claimed: bound + no bond-graph change is a necessary condition only.
  Renamed `intact_bound` (`redox.properties.integrity`, `results/state_integrity.csv`).
- Capacity counted every in-window couple independently and gave ambipolar molecules both
  sides' electrons in either pool. Now: contiguous path from the declared resting state, one
  scorecard row per electrode pool. Counter-ion capacities are scenario conventions.
- Pareto compared only shared objectives, so a candidate with NO clean lambda could dominate
  fully characterized ones (aq_benzylamino dominated aq_benzyloxy and pmdi). Rows missing a
  primary objective are now INCOMPLETE: off the primary front and unable to dominate. The
  front also used sigma_lambda = 0.10 eV while the benchmark gives 0.209 eV; it now uses each
  row's own sigma columns.

**Engineering contracts**
- Energies are protocol-addressed records `calcs/dft/<id>/<state>/sp/<hash>.json`
  (hash = geometry + charge + spin + xc/basis/NLC + solvent + SCF settings); `read_result`
  returns only active-protocol energies and never falls back. Reorg cross-point caches are
  hashed the same way and `_cache_ok` checks the level. `dft.run_batch` flags a result.json
  from another optimization protocol as STALE instead of silently skipping it.
- Missing thermal -> INCOMPLETE (never 0); xTB imaginary modes flagged (`thermal_qc`); thermal
  contribution to each E° reported (`dE_thermal_V`); rotor-cutoff sensitivity in
  `results/thermal_sensitivity.csv`.
- `finalize_after_dft.sh` is fail-fast and now also regenerates dG_disp, its validation, and
  lambda_o.

**Effect on the shortlist.** Anolyte front (both salts): aq_benzyloxy, mophquinone; plus
viologen under TBAPF6. aq_benzylamino and nq_benzyloxy are not dominated but INCOMPLETE (every
lambda couple is QC-flagged: tether conformer jump / unbound gas dianion), so they are not
ranked until a conformer-matched or SMD lambda exists.

**Still open (not fixed by this round):** the E° offset disagreement between OROP and our
anchors (see #23), second-reduction E° validation (1 OROP point), dG_disp validation (n=3),
the gas-phase lambda of dianions (ill-defined when the gas dianion is unbound).

## 23. OROP's +0.5 V E° offset is real; it sits mainly in the Fc reference (see #24)
At the uniform production protocol, the completed OROP subset has a nearly one-sided positive
error (about +0.5 V), while eight in-house events across five families have MAE 0.15 V and
bias +0.13 V. This disagreement must be resolved before claiming absolute E° accuracy or
fitting per-family corrections. Within-charge-class ranking remains useful.

The tempting explanation that OROP uses arbitrary per-solvent experimental reference constants
is **not established**. The repository's `raw_ferrocene-ref-values.txt` calls its entries
"ferrocene simulated values used as reference": they are computed absolute Fc/Fc+ potentials
for each functional/basis/solvent combination (for example B3LYP gives 3.28 V in water,
4.63 V in MeCN, and 4.28 V in DMF). Solvent dependence is physically expected for an absolute
potential. A few combinations are conspicuous outliers (up to 7.04 V), but those values alone
do not show that OROP's experimental values, already tabulated versus Fc, are on inconsistent
scales.

The offset is therefore an **open diagnostic**, not a conclusion about OROP. The next checks
are: read SI Text S2 for the exact thermodynamic cycle and standard-state/reference
conventions; run 3--5 in-house anchors through the OROP geometry path; and run 3--5 OROP
molecules through the full production conformer path. Those crossed calculations separate a
reference/free-energy convention mismatch from geometry/conformer and protocol effects.

## 24. VALIDATION ROUND 3 (2026-10-02) — every decision value re-derived, grounded, and checked
Rule for this round: every reported number must trace to a computation in the repo or a source
that was actually read (DOI + table/page). Results, all reproducible from files listed below:

**Implementation is correct (independent recompute).** `scripts/validation/audit/recompute_axes.py`
re-derives every published axis from raw energy records using only `redox.core.protocol` +
config: **0 discrepancies in 3,391 checks over 9 tables** (`results/validation/audit_recompute.csv`;
now step 12 of the fail-fast finalize chain). Before this round it found stale tables
(ndi_ammonium inputs newer than the tables), the scorecard re-reading 0.1-rounded MW and
half-λ values, an unchecked convergence rule for λ cross points, a truncated Faraday constant
(96485.0; now N_A·e), and a single-label λ QC flag that let `conformer_jump` hide
`anion_unbound` — all fixed.

**Provenance bug in energies (fixed).** 37 anion states (21 pipeline + 16 OROP) optimized before
RI-J became the SCF default (commit fcf05f8, 2026-09-08) were ADOPTED into active-protocol
records that claim density fitting. For the candidates this was 7–12 meV per state
(pmdi, dmophquinone, mophquinone, ndi_ammonium anions). Adoption now requires recorded
density-fitting / tolerance / guess-sweep provenance (`redox.qm.sp._adoptable`,
`redox.core.protocol.record_provenance_ok`); all 85 adopted records were recomputed.
After the fix every candidate state energy equals an independent multi-guess recompute to
≤0.003 meV. Effect on candidates: E1 +4–5 mV, λ_i +7–9 meV, ΔG_disp +0.7–1.0 kJ/mol.

**SCF ground state.** `scripts/validation/reference/scf_level_check.py` re-solved every state
from up to 9 initial guesses (minao/atom/huckel + N±1 orbital-occupation guesses).
Candidates (`results/validation/scf_candidates_check.csv`, current records vs the lowest
solution found): 105 of 106 phases — 63 of 64 state phases and all 42 λ cross points — agree to
≤0.0025 meV. Unchecked: ndi_ammonium red1 in SMD (needs >32 GB GPU; production computed it on
CPU in 12.5 h); its gas phase is on the lowest solution. Reference set (36 states):
two wrong-state production SCFs — **ferrocenium gas phase 2.54 eV high** (E(Fc+)−E(Fc) 8.97 eV
vs NIST evaluated IE 6.71 ± 0.08 eV, webbook.nist.gov CAS 102-54-5; the SMD energy used for E° is
on the lowest solution, so the Fc reference 4.356 V is unaffected) and **OROP system 170
reduced state in SMD, 0.112 eV high**. Open-shell single points now keep the lowest of
minao/atom/huckel (`redox.qm.dft._lowest_scf`). That sweep does NOT reach the Fc+ gas minimum
(atom was still 0.93 eV high; only the occupation guess did): **do not use any gas-phase
quantity of ferrocene** (its IP, ΔG_solv or λ).

**E° against a sourced MeCN benchmark** (`config/benchmark_mecn.py` ←
`data/raw/validation/two_wave_mecn/two_wave_mecn.csv`, 42 rows, every row with DOI + table;
tier A = E°′/E1/2 reported vs Fc directly or with the authors' own calibration; tier B =
half-peak potentials, water-containing electrolyte, or an unverified conversion constant;
excluded = cathodic peaks / irreversible waves / unspecified stereo; 33 new molecules computed
with the production protocol; `results/validation/benchmark_mecn.csv`). Tier A, computed −
experimental:

| couple | n | bias | SD | Spearman |
|---|---|---|---|---|
| E1 (0/−1), all | 24 | +0.26 V | 0.11 V | 0.97 |
| E1, quinones | 18 | +0.31 V | 0.07 V | 0.92 |
| E1, imides | 6 | +0.12 V | 0.10 V | 0.94 |
| E2 (−1/−2) | 8 | −0.07 V | 0.18 V | 0.75 |
| spacing E1−E2 (reference-free) | 9 | +0.25 V | 0.17 V | 0.80 |

Reading: ranking within a family is strong; absolute E1 is systematically too positive
(tight for quinones). The disproportionation over-estimate (ΔG_disp, tier A n=10: MAE
23.4 kJ/mol, all errors positive, RMSE 0.29 eV) comes from E1, not E2 — E2 is close to
experiment (1,2-naphthoquinone is a −0.46 V outlier). A Fc-reference error alone would shift
E1 and E2 equally, so it is not the whole explanation. Nothing is corrected or fitted.

**OROP offset (FINDINGS #23), level-crossing result.** At OROP's own level (B3LYP-D3/6-31G*,
C-PCM ε=37.5; D3 damping not stated by OROP — zero damping assumed) our code gives an absolute
Fc/Fc+ of 4.185 V vs OROP's tabulated 4.662 V, and E vs Fc for 11 OROP systems +0.43 V above
OROP's own B3LYP values (+1/0: +0.38, n=8; 0/−1: +0.57, n=3), while the absolute organic
potentials agree to ~0.07 V (e.g. system 1: 5.615 vs 5.684 V). So most of the OROP gap sits in
the ferrocene reference, not the organics; why OROP's Fc number is higher is not established.
Our Fc gas IE is low vs NIST at both levels (production 6.42 eV incl. xTB ΔZPE, B3LYP 6.19 eV),
consistent with — not proof of — the level-matched Fc reference contributing a positive offset.
`results/validation/{scf_ground_state_check,level_crossing_Efc}.csv`. Refreshed production OROP
benchmark (all records re-verified, n=169): MAE 0.52 V, bias +0.49 V, Spearman 0.92;
within class 0/−1 (n=49) 0.78, +1/0 (n=118) 0.79 (`results/orop_benchmark.csv`).

**Grounded uncertainties now used in the scorecard.** σ_E = per-family RMSE over grounded
points only (benchmark tier A + `grounded=True` anchors; Fc never counted — its 0 residual had
flattered the pooled σ): quinone 0.298 V (was 0.118), imide 0.139 V (was pooled 0.178), viologen
0.088 V (was 0.131). σ_disp 0.288 eV (was 0.183 from 3 uncited points). Anchor audit:
methyl viologen's old values came from an aqueous couple (Bird & Kuhn) — replaced by Cook 2017
(MeCN; absolute values need an unverified Ag/Ag+→Fc constant, spacing 0.420 V grounded);
TEMPO +0.249 V verified (Gerken & Stahl 2015); phenothiazine, anthraquinone and
N-methylpyridinium anchors have no verified source (`grounded=False`).

**Outer-sphere λ_o.** `results/lambda_outer.csv` held only the legacy B3LYP/6-31G(d,p),
non-SCF-checked values although the v2 (def2-SVPD, SCF-checked, both geometries) runs for all 21
candidate couples existed; consolidated now (shifts −8 to −24 meV; aq_benzyloxy and bisviologen gain a
λ_o). `redox.properties.lambda_outer.solvent_constants` imported `config.electrolyte` (never
importable) and silently fell back to hard-coded constants; now loads the config and fails loudly.

**Effect on the shortlist: none.** Pareto-optimal (LiPF6) aq_benzyloxy, mophquinone; (TBAPF6) the
same + viologen — identical before and after every correction above. ndi_ammonium becomes a
complete (dominated) candidate; nq_benzyloxy (no clean λ) leaves the TBAPF6 partial front.
