# TODO

Living task tracker — updated as tasks complete or new ones appear. Strategic phased plan
lives in `docs/PLAN.md`; this is the day-to-day worklist.

Legend: `[ ]` todo · `[~]` in progress · `[x]` done · `[!]` blocked

## Now
- [~] Fix env: reinstall `torch==2.6.0+cu124` with CUDA libs; verify `cuda_avail=True`.
- [ ] Verify fairchem imports and loads UMA (`uma-s-1p2p1`) on GPU.
- [ ] `src/redox/uma.py`: FAIRChem OMol ASE calculator; per-state charge+spin; BFGS
      geometry opt; write relaxed geom + energy to `calcs/uma/<id>/<state>/`. Resumable.
- [ ] P1 smoke test: N-benzylpyridinium ox/red → first adiabatic ΔE; sanity checks
      (no fragmentation, sensible spin).

## Next
- [ ] Parallel UMA runner: fan the 15 states across available GPUs (detect count), skip-if-done.
- [ ] Run UMA gas-phase relaxation + descriptors over the full library.
- [x] Review `docs/DATASETS.md` (research agent) → top DBs identified (below).
- [ ] P2 validation set config: parent redox cores w/ known MeCN vs Fc/Fc⁺ potentials + Fc.

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

## UMA model / driver constraint (confirmed)
- This node's GPU driver = CUDA 12.4 → caps us at `fairchem-core 2.7.1` + `torch 2.6+cu124`
  → only `uma-s-1`, `uma-s-1p1`, `uma-m-1p1` load. `uma-s-1p2`/`1p2p1` fail here
  (`HydraModel` config needs fairchem ≥2.8 → torch ≥2.8 → CUDA ≥12.6).
- **lambda6** has a newer driver (CUDA ≥12.6) and CAN run `uma-s-1p2`/`1p2p1`.
- Plan: develop on this node with `uma-s-1p1`; run production/final on **lambda6** with
  `uma-s-1p2p1` (model is swappable via `--model`). Or ask sysadmin to bump this driver.

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
