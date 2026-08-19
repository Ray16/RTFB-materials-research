# TODO

Living task tracker — updated as tasks complete or new ones appear. Strategic phased plan
lives in `docs/PLAN.md`; this is the day-to-day worklist.

Legend: `[ ]` todo · `[~]` in progress · `[x]` done · `[!]` blocked

## CURRENT STATE (handoff)
- **Env is WORKING**: conda env `redox` = `torch 2.8.0+cu128` + `fairchem-core 2.21.0` +
  pyscf. cu128 (CUDA 12.8) runs on this node's 12.4 driver via CUDA minor-version
  forward-compat (8 GPUs, `cuda=True`). HF logged in (`rayzhu16`), UMA access OK.
- **CORRECTION to earlier notes**: the driver is NOT too old. cu130 (CUDA *13*, major
  mismatch) was the only real blocker; `cu12x` torch runs fine via minor-version compat.
  Production model = **`uma-s-1p2p1`** (newest UMA), enabled via **registration** in
  fairchem 2.21 (`uma.py::ensure_registered` clones the uma-s-1p2 registry entry → 1p2p1
  checkpoint; verified: H2O 1p2p1 vs 1p2 differ ~0.5 kJ/mol, matching lambda6). Direct
  path-load fails (checkpoint's embedded config has newer HydraModel kwargs). DONE:
  requirements (torch2.8/fairchem2.21), setup_env (cu128), uma.py DEFAULT_MODEL=uma-s-1p2p1.
- **Full UMA library was relaxed with uma-s-1p1** (15/15 in `calcs/uma/`). Re-run with the
  new default: remove/rename old `calcs/uma/*/result.json` then `./scripts/run_uma.sh`
  (resumable, so it skips existing results otherwise).

## Now
- [x] Confirm `uma-s-1p2p1` runs (registration verified; H2O cross-check matches lambda6).
- [ ] Re-run UMA library with uma-s-1p2p1, then descriptors on new geoms.
- [ ] Run DFT+SMD batch: `./scripts/run_dft.sh` (all 15 states, CPU fan-out, resumable)
      → then `python -m redox.redox` for the first SOLVATED E° table.
- [ ] Run `python -m redox.descriptors` on final geoms for RMSD (prefers DFT-opt geoms).

## Next
- [x] Parallel UMA runner (`scripts/run_uma.sh`, GPU fan-out) + full library relaxed.
- [x] `docs/DATASETS.md` reviewed → top DBs identified (below); OROP + ReSolvedDB cloned.
- [x] P2 validation set config (`config/validation.py`): 6 parent cores + ferrocene.
- [ ] **DFT+SMD geometry OPTIMIZATION** (not just single points) — `dft.py` currently does
      single points; add SMD opt so λ/RMSD use solvated geoms (per plan).
- [ ] Run validation: parents → pipeline → compare to OROP experimental → MAE → §V gate.
      Replace provisional exp anchors in validation.py with OROP values (physics, not fit).
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
