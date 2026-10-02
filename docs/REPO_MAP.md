# Repo map

A one-screen guide to where things live. See `README.md` for the science, `CLAUDE.md` for
operating rules, `docs/PLAN.md` for the phased plan, `docs/DESIGN_AXES.md` for what each
screening axis means and how far it can be trusted, `docs/DATASETS.md` for external data.

## `src/redox/` — the `redox` package (import as `redox.<subpackage>.<module>`)

`pip install -e .` once (editable install via `pyproject.toml`), **or** prefix commands with
`PYTHONPATH=src`. Every CLI module runs as `python -m redox.<subpackage>.<module>`.

**`core/` — shared helpers and data contracts**
- `common.py` — repo paths, physical constants, `load_config`, `read_manifest`,
  `read_result` (returns ACTIVE-protocol energies only; `sp_status="missing"` otherwise),
  `state_names` (real redox states only — never `reorg/`, `sp/`, `_*`, `*.bak*`),
  `free_energy` (G = E_smd + g_thermal; None if either is missing — never a silent 0).
- `protocol.py` — **the energy protocol** (`ACTIVE_SP`: ωB97M-V/def2-TZVPD for EVERY charge
  state and Fc, SMD MeCN, RI-J, conv 1e-9) and the protocol hash that addresses every energy
  record `calcs/dft/<id>/<state>/sp/<hash>.json` (hash = geometry + charge + spin + level +
  solvent + SCF settings).

**`build/` — structure generation** (SMILES → 3D conformers into `library/`, stereo-checked)
- `structures.py` — shared builder (decorate, conformer ensemble, stereo check, XYZ).
- `candidates.py`, `standalone.py`, `validation.py` — the three library sets (grafted
  candidates, "before functionalization" references, known-E validation cores).

**`qm/` — electronic structure**
- `uma.py` — UMA charge/spin-aware gas pre-opt.
- `dft.py` — DFT+SMD geometry optimization (r2SCAN-D4/def2-SVP(D) in SMD) + single points.
  `run_batch` refuses to silently accept a `result.json` from a different optimization
  protocol (`[STALE]`, needs `--force`). `torsion_scan.py` — rotated seeds for floppy species.
- `sp.py` — compute/adopt the active-protocol energy record for every state
  (`--audit`, `--all --shard n:i`, `--orop` for the benchmark tree).

**`properties/` — per-molecule screening axes**
- `potentials.py` — E° vs a LIVE level-matched Fc/Fc⁺ (same G definition); INCOMPLETE couples
  flagged → `results/redox_potentials.csv`.
- `integrity.py` — bound + intact charged-state screen (`intact_bound`); a necessary
  condition, NOT electrochemical reversibility.
- `stability.py` — disproportionation ΔG (= F·ΔE₁₂).
- `reorg.py` — Nelsen 4-point λ (λ_O, λ_R, sum) with protocol-addressed cross-point caches.
  `nelsen.py` — on-the-fly 4-point for validation scripts.
- `lambda_outer.py` — outer-sphere λ_o (Born, molecular-cavity nonequilibrium PCM; sign-checked,
  point-charge sphere test). `solvated_reorg.py` — SASA radius + Born helper.
- `capacity_and_proxies.py` (capacity, SA score, solubility proxy), `descriptors.py` (RMSD).

**`screening/` — ranking**
- `scorecard.py` — one row per (molecule, electrode pool): contiguous redox path from the
  resting state, explicit λ conventions (λ_het, λ_SE), INCOMPLETE propagation.
- `pareto.py` — σ-aware Pareto front on complete rows only (incomplete rows cannot dominate).

**`validation/` — accuracy checks used by the pipeline**
- `stability.py` — ΔG_disp vs experimental wave spacing → σ_disp.
- `d3tales.py` — D3TaLES dump ingest (shared by validation and candidate mining).

## `scripts/` — entry points (no reusable logic; that belongs in `src/redox/`)
- `check_env.py`, `free_gpus.py` — env check and cluster GPU scan (paths used by `CLAUDE.md`).
- `pipeline/` — the production chain: `run_uma.sh`, `run_dft.sh` (both reserve GPUs via
  `gpu_reserve`), `finalize_after_dft.sh` (fail-fast: tests → SP audit → Fc → E° → λ →
  integrity → stability → λ_o → capacity → scorecard → Pareto → figures), `set_fc_reference.py`,
  `compute_lambda_outer.py` → `consolidate_lambda_outer.py`, `fill_candidates_xlsx.py`.
- `analysis/` — robustness checks on the candidates: `thermal_sensitivity.py` (RRHO cutoff
  spread), `reorg_conformer_matched.py` (tether-rotation-free λ).
- `validation/` — benchmark campaigns, one folder each:
  - `orop/` — `run_orop_benchmark.py` (external E° benchmark, same protocol as production),
    `orop_backfill_thermal.py`.
  - `solvation/` — `validate_solvation.py` (FreeSolv/MNSol ΔG_solv).
  - `reorg_d3tales/` — D3TaLES λ reproduction: set builder, workers at THEIR level
    (`validate_reorg_worker_d3tales.py`) and OUR production level (`…_prod.py`), aggregation +
    unbound-anion QC (`aggregate_d3tales_reorg_prod.py`, `finalize_reorg_qc.py`,
    `reorg_screen_unbound.py`, `reorg_recompute_guarded.py`).
  - `reorg_literature/` — literature λ anchors (oligoacenes, perfluoroacenes):
    `validate_reorg_worker_{generic,inner_acene}.py`, `preflight_reorg_set.py`.
  - `reference/` — `scf_level_check.py`: (1) SCF ground-state check — every state re-solved
    from several initial guesses (`--set reference` = anchors/Fc/OROP sample; `--set
    candidates` = every energy feeding a ranked candidate), (2) level crossing at OROP's
    B3LYP-D3/6-31G* (SMD and C-PCM) → `results/validation/{scf_ground_state_check,
    level_crossing_Efc}.csv`.
  - `audit/` — `recompute_axes.py`: independent from-raw recompute of every published axis
    (imports only `redox.core.protocol` + config) → `results/validation/audit_*.csv`.
- `mining/` — D3TaLES family coverage + low-λ/low-SA candidate mining.
- `fleet/` — multi-node GPU fan-out, every worker through `~/bin/gpu_reserve`
  (`dft_launch.sh` = `redox.qm.dft --all` with one global shard numbering across nodes;
  `sp_fleet_worker.sh` = SP records + λ cross points (`CROSS_IDS`) + OROP SPs;
  `scf_level_launch.sh`; d3level/reorg/OROP workers + launchers; `gpu_probe.sh` +
  `cluster.env` back `free_gpus.py`).
- `polaris/` — ALCF Polaris env build + PBS jobs.
- `plotting/` — **all figures**, one folder per `results/figures/` folder
  (`candidates/`, `validation/`, `reorg/`, `pipeline/`); every `plot_*.py` uses
  `plotting/plot_style.py` `apply_style()`.

## `tests/`
`test_contracts.py` — thermodynamic identities, state ordering, protocol/cache invalidation,
contiguous capacity, Pareto completeness/σ, PCM sign + Born limit.
`PYTHONPATH=src python -m pytest -q tests` (also step T of the finalize chain).

## `config/` — parameters (no logic)
`project.json` + `electrolyte.py` (solvent, window, Fc record), `validation.py` (hand
anchors, each event flagged `grounded`), `benchmark_mecn.py` (sourced MeCN E1/E2/wave-spacing
benchmark from `data/raw/validation/two_wave_mecn/`, tiered A/B/excluded per wave),
`standalone.py`, `starting_candidates.py`, `merrifield_multielectron.py`, `redox_groups.py`,
`scorecard_config.py`.

## Data and outputs
- `data/raw/candidates/` — the input candidate sheets (`Candidates.xlsx`,
  `merrifield_multielectron_smiles.xlsx`). `data/raw/validation/` — external datasets
  (git-ignored) except `two_wave_mecn/`, our hand-curated literature table (tracked; every row
  has DOI + table/page + original reference electrode).
- `library/` — generated structures + `manifest.csv` (resting state = `n_e == 0`).
- `calcs/dft/<id>/<state>/` — `opt.xyz` + `result.json` (geometry, thermal) + `sp/<hash>.json`
  (energies) + `reorg/` cross points. `calcs/spincheck/<id>/<state>_m<k>/` — DFT checks of the
  alternative multiplicity. `calcs/orop/` — benchmark states.
- `results/*.csv` — axis tables (`redox_potentials`, `state_integrity`, `reorganization`,
  `stability_disproportionation`, `lambda_outer`, `scorecard`, `pareto_shortlist`,
  `orop_benchmark`, …), `results/Candidates_computed.xlsx`, `results/figures/{candidates,
  validation,reorg,pipeline}/` (PNG).
- `archive/` — obsolete scripts/modules/outputs kept for provenance (see `archive/README.md`).

## Not in git
bulk `calcs/` outputs, `paper/` (copyrighted PDFs), `data/raw/validation/*/`, `logs/`,
`archive/{tmp,logs}/`, model weights, tokens. See `.gitignore`.
