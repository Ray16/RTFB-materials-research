# Small reproducible CPU example

This example starts from 32 declared design recipes, deduplicates exact molecules,
runs the descriptor stage, replays every output serially, freezes the successor
route stage and verifies that missing exact-route evidence keeps every design on
hold. It also regenerates the checked conditional-grafting figure from its frozen
scientific table and compares the plotted numerical data. It submits no quantum
jobs and requires no CUDA, model weights, raw quantum files or private credentials.

Use Python 3.11 or later with RDKit, NumPy, SciPy, psutil and Matplotlib. Install
the tested direct dependencies with `python -m pip install -r examples/publication_smoke/requirements-cpu.txt`. The exact
tested versions and interpreters are recorded in the execution summary. A separate
plotting interpreter may be supplied when those packages live in two environments.
The full quantum environment remains documented separately in `setup_env.sh` and
`requirements.txt`; this small example does not validate that environment.

From the repository root:

```bash
export PYTHONPATH="$PWD/src"
export OMP_NUM_THREADS=1 OMP_THREAD_LIMIT=1 MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 XTB_THREADS=1
python examples/publication_smoke/run.py --out calcs/publication_smoke
```

The output directory must be new and inside the checkout. Review its
`summary.json`, descriptor and route receipts, and `figures/conditional_grafting.pdf`.
The figure demonstrates reproduction of curated results; it does not recompute
their energies or establish their external accuracy. The molecular subset is an
execution example and carries no synthesis recommendation.

For a separate plotting environment, append `--plot-python /path/to/python`.
