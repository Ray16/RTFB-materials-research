# Scripts

Most scripts here are campaign tools from the project's development: they prepare,
run, collect or review specific calculation campaigns, and many need the raw QM
records and a GPU cluster. The scripts below are the ones behind the JCIM paper.

## Rebuild the paper from committed tables (CPU only, no raw records needed)

| Output | Script | Inputs |
|---|---|---|
| Figure 1 (grafting schematic) | `plotting/manuscript/fig1_grafting.py` | none |
| Figure 2 (pipeline) | `plotting/manuscript/fig2_pipeline.py` | none |
| Figure 3 (QC funnel) | `plotting/manuscript/fig3_qc_funnel.py` | `manuscript/evidence/2026-11-18_cpu_checks/readiness_*.csv` |
| Figure 4a/4b (benchmark parity, residuals) | `plotting/manuscript/fig4_validation_parity.py` | `manuscript/evidence/2026-10-22_benchmark_complete/benchmark.json` |
| Figure 5a/5b (Pareto fronts) | `plotting/manuscript/fig5_pareto.py` | `data/manuscript/full_pool_scored_legacyG_MeCN_TBAPF6.csv` |
| Figure 6a/6b (DMPZ design rule) | `plotting/manuscript/fig6_dmpz_design_rule.py` | same scored table |
| Figure 7 (proposed cells) | `plotting/manuscript/fig7_cell_voltage.py` | same scored table |
| TOC graphic | `plotting/manuscript/fig_toc.py` | same scored table |
| Exact structures of the scored candidates (`data/manuscript/pool_structures.csv`) | `manuscript/pool_structures.py` | `config/*.py` candidate definitions, scored table |
| Supporting Information + `manuscript/si_tables/` | `manuscript/build_si.py` | the records above and `manuscript/evidence/` |
| `manuscript/jcim/manuscript.tex` / `.pdf` | `manuscript/build_latex.py` | `manuscript/manuscript_draft_v0.md`, figure captions |

```bash
for f in scripts/plotting/manuscript/fig*.py; do python "$f"; done
PYTHONPATH=src python scripts/manuscript/build_si.py
python scripts/manuscript/build_latex.py          # needs xelatex + bibtex
```

`validation/check_manuscript_rebuild.py` runs this path plus structure rebuilding
in an isolated export of a committed revision (`--ref HEAD` by default). It
deletes generated artifacts first, compares rebuilt bytes, and preserves logs
and a JSON report. Use `--no-compile` to omit PDF compilation. Dependencies are
listed in `requirements-manuscript.txt`; system Graphviz and DejaVu fonts are
also required. SI software versions come from the committed
`data/manuscript/software_versions.json`, not the rendering environment.

## Recompute the committed records from raw QM records

These read the raw-record archive. Point `REDOX_DATA_ROOT` at the directory that
contains its `calcs/` and `library/`, and put `xtb` on `PATH` (or set `XTB_BIN`).

| Record | Script |
|---|---|
| Scored pool table, verified couple by couple | `validation/score_pool_legacy_tier.py` |
| Benchmark at the pool's xTB thermal tier (`benchmark_legacy_tier.json`) | `validation/benchmark_tier_comparison.py` |
| Seed, outer-sphere, basis and surrogate numbers (`supporting_stats.json`) | `validation/manuscript_supporting_stats.py` |
| Capacity columns of the scored table | `validation/add_capacity_conventions.py` |
| Throughput table (`manuscript/tables/throughput_measured.csv`) | `validation/throughput_report.py` |
| Readiness records (`readiness_*.csv`) | `python -m redox.quality.readiness` |
| Worked-example replays | `validation/check_refactor_oracles.py` |
