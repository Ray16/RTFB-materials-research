# Design axes for grafted multi-electron RFB redoxmers — what we compute, and how well

Target: multi-electron redox-active motifs grafted onto a Merrifield (chloromethyl-polystyrene)
monomer, MeCN / PF₆⁻ electrolyte, Fc/Fc⁺ reference.

Two rankings matter and they do **not** coincide:
- **Importance** — how much the property decides whether the molecule is a viable redoxmer.
- **Reliability** — how well the pipeline computes it *today*, as measured against a benchmark.

Governing principle (FINDINGS #0): energy **differences** where solvation/phase terms cancel are
reliable; anything needing **absolute** charged-species solvation or the **polymer phase** is not.
Many important properties (crossover, swelling, site accessibility, cycle-life chemistry beyond
disproportionation) are **not computable** from single-molecule DFT and are out of scope here.
This file is about getting the **computable** ones right before screening.

Legend: ✅ validated · ⚠️ ranking / relative only · ❌ not reliable or not computed

---

## 1. Accuracy contract for the computable axes (production protocol, re-validated 2026-10-02)

Production protocol: r2SCAN-D4/def2-SVP(D) geometry **optimized in SMD(MeCN)** →
ωB97M-V/def2-TZVPD single points (**one diffuse basis for every charge state and for Fc**,
`redox.core.protocol.ACTIVE_SP`) → G = E_SMD + GFN2-xTB quasi-RRHO → E° vs a **live, level-matched**
Fc/Fc⁺ (same G definition). Every energy is a protocol-addressed record; missing data is
INCOMPLETE, never zero (FINDINGS #22).

| axis | what is computed | benchmark (n) | measured accuracy | status | known limits |
|---|---|---|---|---|---|
| **E°, within-family ranking** | adiabatic ΔG(SMD) per 1e couple | sourced MeCN benchmark, tier A (E1: quinones 18, imides 6); OROP within charge class (0/−1: 49; +1/0: 118) | Spearman 0.92 (quinone E1) / 0.94 (imide E1); OROP 0.78 / 0.79 | ✅ ranking | OROP ranking weaker for 0/−1 couples |
| **E°, absolute** | same, vs level-matched Fc | sourced MeCN benchmark, tier A (35 points, `config/benchmark_mecn.py`) | E1 (0/−1, n=24): bias +0.26 V, SD 0.11 V (quinones +0.31, imides +0.12); all tier A: MAE 0.22 V, bias +0.18 V | ⚠️ systematic offset, family-dependent | OROP +0.5 V offset sits mainly in the Fc reference (level crossing, FINDINGS #24); nothing is corrected; σ_E = per-family RMSE (quinone 0.30, imide 0.14, viologen 0.09 V) |
| **E° of the 2nd reduction** (−1/−2) | same | sourced benchmark tier A (8: benzoquinones, 1,2-naphthoquinone, PMDIs) | bias −0.07 V, MAE 0.12 V, SD 0.18 V | ⚠️ (n=8) | 1,2-naphthoquinone −0.46 V outlier; ion pairing ignored (FINDINGS #19) |
| **Wave spacing / disproportionation ΔG_disp** | F(E₁−E₂) from the same ladder | sourced reference-free spacings, tier A (10: benzoquinones, 1,2-NQ, PMDIs, 2 viologens) | MAE 23.4 kJ/mol, all errors positive, RMSE 0.29 eV, Spearman 0.80 | ⚠️ over-estimated | driven by E1 being too positive (E2 is accurate); stability is OVER-predicted by ~0.25 V |
| **Charged-state integrity** (`intact_bound`) | solvated EA > 0 and no bond-graph change | OROP irreversible cases (CCl₄, CH₂Br₂, anhydride) | removes all known category errors (anion ρ 0.49 → 0.95) | ✅ as a necessary-condition gate | NOT reversibility: misses proton transfer, dimerization, nucleophilic attack |
| **Inner-sphere λ (λ_O, λ_R, sum)** | Nelsen 4-point, gas, SMD geometries | D3TaLES at THEIR level (322) | median \|Δ\| 0.09 eV | ✅ method | production vs D3TaLES: −0.24 eV offset (functional+basis, not error), robust σ 0.21 eV (330) |
| **λ QC** | conformer-jump / unbound-anion / negative-half flags | thiosuccinimide outliers (16/442) | all caught | ✅ | flagged couples carry no λ → candidate INCOMPLETE on kinetics |
| **Outer-sphere λ_o** | 1-body nonequilibrium PCM (molecular cavity) + Born | point charge in a spherical cavity | PCM = Born to <1e-5 | ⚠️ ranking | continuum ceiling ~0.3 eV vs explicit solvent; not yet checked against a second PCM code |
| **λ_het, λ_SE(d)** | (λ_O+λ_R)/2 + λ_o,1 ; λ_O+λ_R + 2λ_o,1(1−a/d) | — | derived | ⚠️ annotation | inherits λ_o limits; polymer site distance d not known |
| **Capacity** | n F / mass on the contiguous path from the resting state, per electrode | — | bookkeeping | ✅ | counter-ion variants are scenario conventions |
| **Spin ground state** | UMA multiplicity scan; DFT check when the gap is small | aq_benzyloxy dianion | UMA gap 0.19 eV vs DFT 0.81 eV | ⚠️ | UMA underestimates gaps; DFT-check every gap < 0.5 eV |
| **SCF solution** | open shells: lowest of minao/atom/huckel | up-to-9-guess recompute of every candidate energy (64 state phases + 42 cross points) | 105/106 on the lowest solution found (≤0.0025 meV); ndi_ammonium red1 SMD unchecked (CPU-only size) | ✅ | ferrocenium GAS SCF is 2.5 eV high (E° unaffected; never use Fc gas quantities) |
| **Thermal (RRHO)** | GFN2-xTB qRRHO at the DFT-SMD geometry | rotor cutoff 25/50/100 cm⁻¹ | see `results/thermal_sensitivity.csv` | ⚠️ | some states have imaginary xTB modes (flagged `thermal_qc`) |
| **ΔG_solv (neutral)** | SMD | FreeSolv + MNSol (31 neutrals) | MAE 0.79 kcal/mol | ✅ | relative only for solubility; does not transfer monomer → polymer |
| **SA score** | Ertl | — | proxy | ⚠️ filter | does not measure graftability |

Figures: `results/figures/validation/orop_benchmark_parity.png` (E°),
`validation.png` (in-house anchors), `reorg/reorg_prod_parity_3panel.png` and
`reorg/d3level_hole_electron_parity.png` (λ), `solvation_parity.png` (ΔG_solv).

## 2. What must be true before screening at scale (ordered)
1. **Resolve the E° offset** between OROP and the in-house anchors (§1). Leading hypotheses:
   a reference/free-energy convention mismatch versus a protocol/geometry difference (OROP
   systems are re-optimized from supplied geometries without xtb-ALPB pre-opt or
   multi-conformer seeding). OROP's `raw_ferrocene-ref-values.txt` contains **simulated absolute
   Fc potentials**, separately for each solvent and electronic-structure method; their
   solvent dependence is expected and is not evidence that the experimental Fc scale is
   inconsistent. The large outliers across some method/solvent combinations are a reason to
   audit SI Text S2, not an explanation by themselves (FINDINGS #23).
   Test: run 3–5 anchors through the OROP path and 3–5 OROP molecules through `run_batch`.
2. **Per-family E° calibration** only after (1), with ≥3 experimental anchors per family
   (imides currently have none).
3. **Validate second reductions** (E₂) and grow the ΔG_disp benchmark to ≥10 molecules.
4. **Make λ computable for every candidate couple**: conformer-matched core λ for tether jumps;
   SMD (not gas) λ for dianions whose gas state is unbound.
5. DFT spin checks for all small UMA gaps; thermal sensitivity below the E° σ.

---

## 3. Tiers (importance) — unchanged in substance

### Tier 1 — decisive
| axis | why it matters | status |
|---|---|---|
| **Redox potential E°** | energy density ∝ n·ΔV·C | ⚠️ ranking; absolute unresolved (§1–2) |
| **Multi-electron wave spacing / window fit** | all n electrons usable at one voltage inside the window | ⚠️ ΔG_disp over-predicted by ~24 kJ/mol (n=10); window gate on the contiguous path |
| **Charged-state stability** | #1 failure mode of organic RFBs | ✅ disproportionation · ❌ dimerization/pimerization (continuum fails for +2 desolvation) · ❌ general decomposition (needs a per-family reaction library) |
| **Capacity: n per total mass** | grafting dilutes gravimetric capacity | ✅ bookkeeping (monomer, repeat unit, counter-ion scenarios) |

### Tier 2 — architecture-dependent
| axis | status |
|---|---|
| **Integrity gate** (`intact_bound`) | ✅ necessary-condition pre-filter, not a ranking axis |
| **Solubility** (only if the polymer flows) | ⚠️ relative ΔG_solv proxy; ❌ monomer → polymer transfer |
| **Solid architecture: accessibility, swelling, conductivity** | ❌ not computed; hopping handle = λ_SE(d) + transfer integrals |
| **Crossover** | ❌ not computable here; addressed by grafting (immobilization) |

### Tier 3 — secondary
| axis | status |
|---|---|
| **Inner-sphere λ** | ✅ method-validated; rarely limiting |
| **Outer-sphere λ_o, λ_het, λ_SE** | ⚠️ ranking annotation |
| **SA score** | ⚠️ filter |
| **Grafting compatibility** | modeled via the 4-methylbenzyl proxy; parent-vs-grafted shifts via `_sa` references |

---

## 4. Where each axis lives in the output

`results/scorecard.csv` — **one row per (molecule, electrode pool)**:

| axis | column(s) |
|---|---|
| status | `status` (candidate / REJECTED / INCOMPLETE), `reason`, `path`, `path_stop`, `resting_state` |
| E° | `E_V` (pool mean over the path), `E_anolyte_V` / `E_catholyte_V`, `sigma_E_V` — raw, uncalibrated |
| capacity | `n_accessible`, `specific_capacity_mAh_g`, `capacity_repeat_mAh_g`, `capacity_maxload_{Li,TBA}_mAh_g` |
| stability | `dG_disp_kJmol`, `sigma_disp_eV`, `disp_applicable` (False for a 1e path) |
| inner λ | `lambda_i_eV` (λ_O+λ_R), `lambda_i_ox_eV`, `lambda_i_red_eV`, `sigma_lambda_eV`, `lambda_qc`, `lambda_flags`, `n_lambda_flagged`, `n_lambda_missing` |
| outer / combined λ | `lambda_o_pcm_eV`, `lambda_o_method`, `lambda_het_eV`, `sigma_lambda_het_eV`, `lambda_se_contact_eV` |
| QC | `all_intact_bound`, `thermal_qc`, `min_spin_gap_eV`, `spin_gap_source`, `spin_confidence` |
| proxies | `SA_score`, `dGsolv_proxy_eV` |

Per-couple detail: `redox_potentials.csv` (`status`, `thermal_qc`, `dE_thermal_V`),
`state_integrity.csv`, `reorganization.csv`, `stability_disproportionation.csv`,
`lambda_outer.csv` (with `method`, `scf_checked`), `multielectron_stability.csv`.

## 5. Pareto definition (`src/redox/screening/pareto.py`)
- Pools: anolyte / catholyte (a molecule enters a pool only through its own contiguous path).
- Objectives (σ-aware, higher = better): voltage, capacity, stability (ΔG_disp; not applicable
  for a 1e path), kinetics (−λ_i, QC-clean only). σ = **each row's own** scorecard σ.
- **Completeness**: a row missing any applicable objective is INCOMPLETE — never on the primary
  front and never able to dominate; listed as `partial_front` if nothing complete dominates it.
- Capacity is evaluated under each salt scenario (`CAPACITY_SCENARIOS`: LiPF₆, TBAPF₆).
- Annotations only: λ_o, λ_het, λ_SE, SA, ΔG_solv. No scalarized figure of merit.
- Current anolyte front (both salts): **aq_benzyloxy, mophquinone**; **viologen** also under
  TBAPF₆. Not dominated but INCOMPLETE (no QC-clean λ): **aq_benzylamino, nq_benzyloxy**.

Figures (`results/figures/candidates/`): `pareto_parallel_coords.png` (incomplete rows dotted),
`pareto_scatter_matrix_{LiPF6,TBAPF6}.png`, `pareto_dominance_heatmap.png`,
`merrifield_summary.png`, `candidate_E_redox_landscape.png`.
