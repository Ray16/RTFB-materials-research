# Plan — next steps

Phased plan for Stage 1 (decorate Cl site → stability + redox potential). Validation is
woven in, not deferred (§V). Status: **P0 mostly done**, env finishing.

## P0 — Scaffold & conformers  ✅ (env pending)
- [x] Repo + GitHub, folder structure, CLAUDE.md/README/MODELING.
- [x] `config/redox_groups.py` (6 groups, 15 states) + `config/electrolyte.py`.
- [x] `build.py`: decorate Cl site, charge-aware conformer ensemble → `library/`.
- [ ] `redox` env finishes (torch+fairchem) → `python scripts/check_env.py` green.

## P1 — UMA online (gas-phase engine)
- [ ] `huggingface-cli login`; accept `facebook/UMA` license. **← needs user token.**
- [ ] `src/redox/uma.py`: FAIRChem ASE calculator (OMol task), sets **charge+spin** per
      state, BFGS geometry opt, writes relaxed geom + energy → `calcs/uma/<id>/<state>/`.
- [ ] Smoke test: N-benzylpyridinium ox/red → first adiabatic ΔE; confirm charge/spin
      handling and no fragmentation.
- [ ] Parallel runner: fan states across the 4 GPUs (resumable, skip-if-done).

## P2 — Validation set (do BEFORE trusting decorated monomers)  ← §V
- [ ] `config/validation.py`: parent redox cores with **known experimental E° in MeCN vs
      Fc/Fc⁺** (methyl viologen, TEMPO, phenothiazine, 9,10-anthraquinone, N-methyl-
      pyridinium) + ferrocene (internal reference).
- [ ] Run full pipeline on them; build calibration + uncertainty (§V).

## P3 — Solvated optimization (final geometries)
- [ ] Add `xtb` (+ALPB MeCN) and `gpu4pyscf` to the env.
- [ ] xtb-ALPB pre-opt for **charged** states (warm-started from UMA/ETKDG) — kills
      gas-phase salt-bridge artifacts.
- [ ] DFT + SMD(MeCN) geometry opt + frequencies (gpu4pyscf on V100s) → final geoms +
      `G_solv` → `calcs/dft/`. One job per GPU.

## P4 — Descriptors & redox potentials
- [ ] `src/redox/redox.py`: assemble `G_solv`, E° per 1e event, reference to Fc/Fc⁺.
- [ ] `src/redox/descriptors.py`: λ (4-point Marcus), RMSD + structural change, HOMO/LUMO,
      IP/EA, spin-density localization, partial charges, connectivity/stability check.
- [ ] Aggregate → `results/` tables + electronic-property distribution plots.

## P5 — Analysis & report
- [ ] Rank groups by E°, stability, # accessible electrons; flag multi-electron winners.
- [ ] Compare against D3TaLES; report predicted E° **with error bars** (§V).

---

## §V — Validation strategy (we are using approximations, so anchor them)

**External (against reality) — the strongest checks:**
1. **Experimental parents.** Our decorated monomers are benzyl-substituted versions of
   well-characterized cores. Compute the bare cores and compare E° to tabulated **MeCN vs
   Fc/Fc⁺** values (methyl viologen ≈ −0.45/−0.89 V vs SCE; phenothiazine •+ ≈ +0.25;
   TEMPO oxoammonium ≈ +0.20…+0.30; anthraquinone 1st red ≈ −1.2…−1.4). Substituent shift
   from core→benzyl should be small and in the Hammett-expected direction.
2. **D3TaLES.** Pull overlapping/analogous molecules (viologens, phenothiazines, TEMPO,
   quinones, pyridiniums) and compare our E° to their reported DFT + experimental values.
   Optionally match their exact protocol for a direct number-to-number check.

**Calibration & uncertainty:**
3. Fit `E_exp = a·E_calc + b` over the validation set → MAE / systematic offset (redox-DFT
   typically ~0.1–0.2 V); apply the fit to decorated monomers and report residual as the
   error bar.
4. **Method sensitivity:** vary functional/basis on 2–3 systems → spread = method uncertainty.

**Internal (self-consistency):**
5. **UMA vs DFT cross-check:** same subset both ways → geometry RMSD + ΔE correlation;
   validates using UMA geometries and gives the UMA→DFT offset.
6. **Internal Fc/Fc⁺ reference** computed identically cancels systematic error.
7. **Physical sanity:** true minima (no imaginary freqs), ⟨S²⟩ near ideal (low spin
   contamination), connectivity preserved (no fragmentation), Rg (no gas-phase collapse),
   ordering matches the known electrochemical series (n-type ≪ p-type).
8. **Thermodynamic consistency:** sequential 1e potentials vs 2e; disproportionation.
