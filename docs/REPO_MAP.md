# Repo map

A one-screen guide to where things live. See `README.md` for the science, `CLAUDE.md` for
operating rules, `docs/PLAN.md` for the phased plan, `docs/DATASETS.md` for external data.

## `src/redox/` — reusable pipeline code (import as `redox.<module>`, needs `PYTHONPATH=src`)

**Structure generation**
- `build.py`, `build_validation.py`, `build_candidates.py`, `build_standalone.py` — SMILES → 3D
  conformers into `library/` (stereo-checked).
- `uma.py` — UMA (MLIP) charge/spin-aware gas pre-optimization.

**Electronic structure & properties**
- `dft.py` — DFT+SMD driver: in-solvent geometry optimization + gas/SMD single points
  (r2SCAN-D4/def2-SVP(D) // ωB97M-V/def2-TZVP(D), RI-J on). `run_batch` (multi-conformer,
  optional `--torsion-scan`), `dft_smd`. **The production compute path.**
- `reorg.py` — inner-sphere λ_i (Nelsen 4-point) over the pipeline's cached SMD energies;
  artifact guard (`flag`). `nelsen.py` — the same 4-point but for on-the-fly GAS recompute
  (validation/diagnostics), shared by the D3TaLES validation scripts.
- `solvated_reorg.py` — outer-sphere λ_o (Born/Marcus) + r_hyd (SASA).
- `torsion_scan.py` — detect the primary ring-ring torsion; rotated seeds for `--torsion-scan`.
- `redox.py` — redox potentials (Fc/Fc⁺-referenced). `d3tales_ingest.py` — D3TaLES ingest.
- `descriptors.py`, `capacity_and_proxies.py`, `stability.py`, `reversibility.py`,
  `ionpair.py`, `dimerize.py`, `microsolvate.py` — descriptors & stability proxies.
- `scorecard.py`, `pareto.py` — unified scorecard + Pareto shortlist.

## `scripts/` — runnable entry points
- `scripts/plotting/` — **all figures.** Every `plot_*.py` imports the co-located `plot_style.py`
  (`apply_style()`) — the one canonical publication style (18 pt, no overlaps, top-journal theme,
  Okabe-Ito). Run e.g. `python scripts/plotting/plot_results.py`.
- `scripts/diagnostics/` — one-off diagnostics kept for provenance: `diag_gas_anion_reorg.py`,
  `diag_smiles_reorg.py`, `scan_dihedral_reorg.py`, `aggregate_scan_reorg.py`.
- Top-level `scripts/` — the interconnected pipeline/launch/finalize chain (kept flat because
  shell scripts and `CLAUDE.md` reference these by path):
  - launch/env: `run_*.sh`, `free_gpus.py`, `check_env.py`.
  - D3TaLES reorg validation: `build_d3tales_reorg_set.py` → `validate_reorg_worker.py`
    (resumable; `_launch_node_workers.sh` / `launch_reorg_validation.sh` for cluster fan-out) →
    `aggregate_d3tales_reorg.py`; figure via `scripts/plotting/plot_d3tales_reorg_compare.py`.
  - finalize chain: `wait_and_finalize.sh` → `finalize_after_dft.sh` → `set_fc_reference.py` +
    `scripts/plotting/plot_results.py`; plus `add_thermal.py`, `probe_*`, `sample_dimer.py`, etc.

## `config/` — parameters (no logic)
`electrolyte.py` + `project.json` (solvent/referencing), `validation.py`, `standalone.py`,
`starting_candidates.py`, `scorecard_config.py`.

## `results/` — small, committed artifacts (numbers + figures; bulk calcs are git-ignored)
- `*.csv` — property tables (`reorganization.csv`, `scorecard.csv`, `redox_potentials.csv`, …).
- `figures/` — **all figures, in topical subfolders**: `candidates/`, `validation/`, `reorg/`,
  `pipeline/`. Slides find them via `\graphicspath`.
- `d3tales_reorg_validation/` — validation DATA (`comparison.csv`, `molecules.csv`); the
  per-molecule `calc/` JSONs are regenerable and git-ignored.
- `dihedral_scan/` — the two biaryl-quinone torsion scans.

## Not in git
`calcs/` (bulk DFT/UMA outputs), `paper/` (copyrighted PDFs), `data/raw/validation/*/` (re-clone
via `docs/DATASETS.md`), model weights, tokens. See `.gitignore`.
