# TODO

Living task tracker — updated as tasks complete or new ones appear. Strategic phased plan
lives in `docs/PLAN.md`; this is the day-to-day worklist.

Legend: `[ ]` todo · `[~]` in progress · `[x]` done · `[!]` blocked

## CURRENT STATE (handoff) — updated 2026-08-21

### Latest session (thermal + ion-pair verdict + OROP benchmark)
- **THERMAL free-energy corrections: DONE, whole table.** Uses **GFN2-xTB RRHO** (298.15 K)
  via `dft._thermal_correction` (xtb installed; in setup_env.sh + check_env.py). We do NOT use
  the DFT Hessian: **gpu4pyscf's UKS analytic Hessian is numerically BROKEN** in this version
  (inflates open-shell frequencies ~2x → corrupt ZPE, e.g. 14.7 vs correct 6.4 eV). xtb thermal
  is fast, correct for open+closed shell, method-transferable. G = E_smd + g_thermal everywhere.
  Effect: viologen bare-ion MAE 0.429 → **0.372**.
- **ION-PAIR (released-counterion PF6-) scheme: REFUTED** in implicit solvent. Even with a
  correct singlet mv_ip2, waves/spacing are catastrophically wrong (spacing ~+11 V vs exp -0.43):
  forming a neutral ion pair from two well-solvated ions costs a huge (~+6 eV), poorly-modeled
  desolvation free energy that does NOT cancel between waves. See `src/redox/ionpair.py` CONCLUSION.
  BEST viologen = 'bare ions + thermal' (MAE 0.372). Do not resurrect continuum ion-pairing.
- **OROP experimental benchmark: DONE (11/14).** `scripts/run_orop_benchmark.py` → 
  `results/orop_benchmark.csv`. Our **MAE 0.869 V** (≈ OROP raw implicit-DFT 0.821 on same
  systems), **Spearman 0.89** (ranking reliable), **+0.77 V systematic bias** dominated by
  charged species (radical anions +1.24, dication +1.07). De-biased MAE 0.48. Confirms the
  ~0.5-0.9 V floor of pure-physics implicit DFT; the field closes the gap with a charge-dependent
  ML/linear correction (OROP's "implicit+ML") — NOT more physics.
- **NEXT (recommended):** charge-class linear calibration E_exp ≈ a·E_calc + b fit on OROP with a
  held-out split → expected MAE ~0.3-0.5. Also: expand OROP subset (30-60, all charge classes)
  and fix 3 SMD gradient-convergence failures (sys 107/109/115: "Nuclear gradients not converged").

### Earlier state

- **Env is WORKING**: conda env `redox` = `torch 2.8.0+cu128` + `fairchem-core 2.21.0` +
  pyscf 2.14 + **geometric 1.1.1** (geomeTRIC, in-solvent DFT optimizer). cu128 (CUDA 12.8)
  runs on this node's 12.4 driver via CUDA minor-version forward-compat (8 GPUs).
- Production model = **`uma-s-1p2p1`** (newest UMA), enabled via **registration** in
  fairchem 2.21 (`uma.py::ensure_registered` clones uma-s-1p2 registry entry → 1p2p1).
- **UMA spin/conformer scan (uma-s-1p2p1) RUNNING** on GPU1 → `calcs/uma/<id>/<state>/`
  (log `calcs/uma/scan.log`). Scans conformers×multiplicities, picks lowest-E spin state,
  records S-T gap. Verified states so far: pyridinium ox(m=1,ST 3.51eV)/red(m=2),
  cyanopyridinium ox/red(m=2), viologen in progress.
- **DFT now OPTIMIZES IN SOLVENT** (SMD gradients + geomeTRIC), not single-point — matches
  the spec ("geometries optimized in solvent"). Writes `calcs/dft/<id>/<state>/opt.xyz`.
  DFT reads the UMA-chosen multiplicity (`_uma_mult`), not the manifest hint.
- **Referencing upgraded to level-matched ferrocene** (physics not fitting): redox.py
  prefers `FC_ABS_COMPUTED_V` (our own Fc/Fc+ at our level) over the thermodynamic SHE+Fc
  constant. OROP MeCN Fc refs stored in electrolyte.py for cross-check.
- **§V gate re-benchmarked against OROP 313 (MeCN)**: raw implicit-DFT MAE≈0.50 V (bias
  +0.31 V) — sub-0.15 V absolute from pure physics is not realistic. New gate judges
  absolute (≤0.25 good after Fc referencing) AND ranking (Spearman ρ≥0.9) separately.

## Now
- [x] Confirm `uma-s-1p2p1` runs (registration verified; H2O cross-check matches lambda6).
- [x] **DFT+SMD geometry OPTIMIZATION** wired into dft.py (geomeTRIC + SMD gradients).
- [x] UMA spin/conformer scan with uma-s-1p2p1 (all states scanned).
- [x] Full DFT+SMD-opt batch (18 states) → `redox.redox` solvated E° table
      (`results/redox_potentials.csv`).
- [x] Descriptors on DFT-opt geoms → `results/structure_descriptors.csv` (RMSD).
- [x] Validation gate: **MAE 0.18 V, RMSE 0.25 V, mean signed −0.09 V** (n=8, 4/8 within ±0.15 V).
- [x] **Publication figure set** (300 dpi, one panel/figure, 18 pt, no overlap):
      `pipeline.png`, `validation.png`, `redox_landscape.png`, `structure_change.png`;
      stale multi-panel figures pruned.
- [x] **2D structure galleries** (`scripts/plot_molecules.py`, driven off manifest.csv):
      `molecules_candidates.png` (6 monomers) + `molecules_validation.png` (6 cores).
- [x] **Pipeline diagram rebuilt in Graphviz** (`scripts/plot_pipeline_graphviz.py`,
      auto-layout, HTML records) → refreshed `pipeline.png`; graphviz registered in
      requirements.txt + check_env.py (system `dot` verified present).

## Next
- [x] Parallel UMA runner (`scripts/run_uma.sh`, GPU fan-out) + full library relaxed.
- [x] `docs/DATASETS.md` reviewed → top DBs identified (below); OROP + ReSolvedDB cloned.
- [x] P2 validation set config (`config/validation.py`): 6 parent cores + ferrocene.
- [x] **DFT+SMD geometry OPTIMIZATION** — done (geomeTRIC + SMD gradients in dft.py).
- [ ] Run validation: parents → pipeline → compare to OROP experimental → MAE → §V gate.
      Replace provisional exp anchors in validation.py with OROP values (physics, not fit).
- [ ] Compute our own ferrocene Fc/Fc+ reference at our level → set FC_ABS_COMPUTED_V.
- [ ] Ferrocene needs a metallocene geometry (RDKit can't embed) — special-case it.

### Validation / HTS data sources (from docs/DATASETS.md)
- [ ] **OROP/ROP313** — clone; experimental MeCN ox/red potentials → §V calibration anchor
      (`E_exp = a·E_calc + b`). Verify license before redistributing.
- [ ] **D3TaLES** — download no-login CC-BY MDF/Globus dump (DOI 10.18126/v5sj-6q93);
      filter by motif SMARTS; pull MeCN oxidation/reduction potentials + λ + solubility.
- [ ] **ReSolvedDB** — clone; 19,785 computed MeCN reduction potentials → computed cross-check.
- [ ] **OMol25** — HTS candidate pool (HF-gated, aligns with UMA); derive potentials ourselves.

## Soon
- [ ] Add `xtb` (+ALPB MeCN) and `gpu4pyscf` to env (update requirements.txt + check_env.py).
- [ ] xtb-ALPB solvated pre-opt for charged states (warm-start from UMA/ETKDG).
- [ ] DFT+SMD(MeCN) geometry opt + frequencies → `G_solv`.
- [ ] `src/redox/redox.py`: E° per 1e event, referenced to Fc/Fc⁺.
- [ ] `src/redox/descriptors.py`: λ (4-point), RMSD, structural change, HOMO/LUMO, spin
      density, connectivity/stability.
- [ ] Calibration vs experiment/D3TaLES; report E° with error bars.

## Decisions pending
- [ ] Level of theory for DFT+SMD (match D3TaLES vs OMol-aligned range-separated hybrid).
- [ ] Thermal corrections: full harmonic freq vs xtb vs electronic-only first pass.
- [ ] UMA model size for production: `uma-s-1p2p1` vs `uma-m-1p1` (decide after P2 benchmark).

## Done
- [x] Repo scaffold + GitHub (`Ray16/RTFB-materials-research`).
- [x] `config/redox_groups.py` (6 groups, 15 states) + `config/electrolyte.py` (MeCN/SMD/Fc).
- [x] `build.py`: decorate Cl site, charge-aware conformer ensemble → `library/` + manifest.
- [x] Design docs: `docs/MODELING.md`, `docs/PLAN.md` (with validation strategy).
- [x] Reproducible env: `setup_env.sh` (GPU/CPU auto-detect), `requirements.txt`.
- [x] HF auth + UMA gated access confirmed.
