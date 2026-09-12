# scripts/archive — completed one-off experiments

These scripts drove experiments that are **finished** and whose conclusions now live in
`FINDINGS.md`. They are kept for provenance/reproducibility, not for the live pipeline. None
is referenced by the current pipeline chain (see `docs/REPO_MAP.md`).

| script | what it did | outcome (see FINDINGS.md) |
|--------|-------------|---------------------------|
| `run_orop_benchmark.py`      | OROP 313-system MeCN redox benchmark | accuracy floor ~0.5 V (#1) |
| `add_thermal.py`, `run_thermal_batch.py` | batch GFN2-xTB RRHO thermal corrections | now folded into `dft.py` `do_freq` (#4) |
| `run_viologen_fix.py`, `finalize_viologen_fix.py` | viologen re-optimization pass | superseded by `--torsion-scan` in `dft.py` |
| `sample_dimer.py`, `finalize_dimer_sampling.py` | viologen π-dimer sampling | dimerization unreliable in continuum (#9) |
| `finalize_microsolv.py`      | explicit 4-MeCN microsolvation finalize | no accuracy gain, dropped (#6) |
| `probe_thermal.py`           | one-off DFT vs xTB Hessian probe | GPU UKS Hessian broken (#3) |

To re-run one, invoke it from the repo root with the `redox` env active (same as any pipeline
script). They still `import redox.*`, so no path changes are needed.
