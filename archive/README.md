# archive/ — obsolete or completed material, kept for provenance

Nothing here is part of the live pipeline, and nothing in `src/`, `scripts/` or `config/`
imports from it. Conclusions from these experiments live in `FINDINGS.md`. Files keep their
original relative paths under this directory (e.g. `archive/scripts/diagnostics/…` was
`scripts/diagnostics/…`).

| path | what it is | why archived |
|---|---|---|
| `scripts/` | one-off experiment drivers (thermal batch, viologen fix, dimer sampling, microsolvation finalize, thermal probe), batch launchers for finished campaigns (`run_*_batch2.sh`, `run_dft_valcore.sh`, `run_uma_bygroup.sh`, `_launch_node_workers.sh`, `launch_reorg_validation.sh`), env probes (`probe_dft.py`, `probe_uma.py`), and the superseded reorg validation scripts (`validate_reorg_worker.py`, `aggregate_d3tales_reorg.py`, `diagnostics/*`, `plotting/plot_d3tales_reorg_compare.py`, `plotting/plot_starting_candidates.py`) | finished, or superseded by the scripts under `scripts/reorg/` and `scripts/pipeline/` |
| `scripts/fleet/{gpu_guard,gpu_enforce_host,gpu_watchdog}.sh` | per-repo GPU lock/enforcer/watchdog | superseded by the system-wide `~/bin/gpu_reserve` gate (watchdog was already off) |
| `scripts/fleet/d3level_{rolling,refresh}.sh` | hourly monitor/re-aggregator for the D3TaLES-level sweep | sweep finished (468/468); re-aggregate with `scripts/validation/reorg_d3tales/aggregate_d3tales_reorg_prod.py` + `finalize_reorg_qc.py` |
| `scripts/pipeline/wait_and_finalize.sh` | unattended finisher that auto-committed/pushed | auto-commit conflicts with the review-before-commit workflow; run `finalize_after_dft.sh` directly |
| `scripts/plotting/plot_pipeline_matplotlib.py` | hand-laid matplotlib workflow diagram | superseded by the Graphviz version, now `scripts/plotting/pipeline/plot_pipeline.py` |
| `src/redox/dimerize.py` | π-dimer builder | only used by the archived dimer-sampling experiment |
| `results/viologen_ionpair.csv` | ion-pair E° table | output of the archived `ionpair.py` (NEGATIVE result, FINDINGS #5) |
| `src/redox/ionpair.py` | released-counterion ion-pair redox scheme | NEGATIVE result (FINDINGS #5): continuum ion pairing fails |
| `src/redox/microsolvate.py` | explicit 4-MeCN cluster-continuum builder | NEGATIVE result (FINDINGS #6): no accuracy gain |
| `results/reversibility.csv` | old gate output with a `reversible` verdict | renamed to `results/state_integrity.csv` (`intact_bound`), which does not over-claim electrochemical reversibility |
| `results/d3tales_reorg_validation/`, `results/dihedral_scan/`, `results/figures/…` | outputs of the superseded B3LYP/6-31G* reorg comparison and torsion scans | superseded by the D3TaLES-level and production-level comparisons |
| `docs/TODO_2026-08.md` | the August handoff worklist | replaced by a fresh `TODO.md` |
| `slides_build/` | LaTeX build artifacts (`.aux/.nav/.out/.snm/.toc/.log`) | regenerable; now git-ignored |
| `tmp/`, `logs/` | scratch scripts and logs of finished runs | git-ignored; kept on disk only |

Archived Python uses the OLD flat module names (`redox.dft`, `redox.common`, …). The live package
was reorganized into subpackages on 2026-10-02 — map old → new as
`common, protocol → redox.core.*` · `build → redox.build.structures`,
`build_{candidates,standalone,validation} → redox.build.{candidates,standalone,validation}` ·
`uma, dft, sp, torsion_scan → redox.qm.*` · `redox → redox.properties.potentials`,
`integrity, stability, reorg, nelsen, lambda_outer, solvated_reorg, capacity_and_proxies,
descriptors → redox.properties.*` · `scorecard, pareto → redox.screening.*` ·
`validate_stability → redox.validation.stability`, `d3tales_ingest → redox.validation.d3tales`.
To rerun an archived script, update its imports with that map; the archived `src/redox`
modules (`ionpair`, `microsolvate`, `dimerize`) need `PYTHONPATH=archive/src:src`.
