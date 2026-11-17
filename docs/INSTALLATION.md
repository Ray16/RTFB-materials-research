# Installation and execution environments

## Read a cached result without installing the package

From the checkout root:

```bash
python3 -I -S examples/cached_thermochemistry/reproduce.py
```

Python ≥3.10 and the standard library suffice. The supplied files are the inputs;
the result is JSON on stdout, optionally saved with `--out PATH`. This is the
smallest verified example and does not launch a quantum calculation.

## Inspect the package

Use the checkout in place. Its Python package discovers `config/`, `library/`
and other project directories relative to the source tree. A standalone wheel
is not currently a self-contained scientific installation.

```bash
python -m pip install --no-deps --no-build-isolation -e .
python -m redox --help
python -m redox workflows --json
```

This needs an existing pip/setuptools environment. It does not install scientific
packages. For a temporary source-path invocation without installing anything:

```bash
PYTHONPATH=src python -m redox workflows
```

## Scientific calculations

`requirements.txt` records the Python scientific stack. `setup_env.sh` creates a
Conda environment and installs CUDA-matched PyTorch/GPU4PySCF separately. Read
the script before running it: it creates an environment and downloads large
packages, and it is not an automatic step in the quick start.

```bash
CUDA_TAG=cu128 bash setup_env.sh redox
conda activate redox
python -m pip install --no-deps --no-build-isolation -e .
```

The full setup has not been recreated from scratch in this repository-polish
pass. The working runtime is Python 3.11 with the pinned scientific versions;
platform availability, driver compatibility and fresh installation require their
own verification. CUDA package resolution can change even with top-level pins.
`CUDA_TAG=cpu` selects a CPU PyTorch bootstrap and skips GPU4PySCF; this does not
qualify GPU/Hessian campaigns on a CPU-only host.

| Work | Dependencies and constraints |
| --- | --- |
| CLI help and workflow catalogue | Standard library |
| Cached PMDI arithmetic example | Standard library, supplied cached data |
| Structure construction | RDKit and the dependencies imported by the chosen builder; `requirements.txt` records the scientific versions |
| UMA warm starts | PyTorch, fairchem, model access/weights and qualified backend |
| DFT/native thermochemistry | PySCF, GPU4PySCF for GPU work, geomeTRIC, xTB 6.7.1 integrator and declared runtime profile |
| Manuscript builds | pdfLaTeX, BibTeX, pypdf and the self-contained manuscript inputs |
| Full test suite | pytest plus scientific dependencies exercised by the selected tests |

UMA model access and third-party data may require separate authorization. Never
place credentials in source/configuration files. A missing reservation service
must not be bypassed by setting CUDA visibility manually. The shared resource
layer supports declared gate discovery; consult `docs/WORKFLOWS.md` before use on
another cluster.

## Inputs and outputs

Candidate configuration files are executable, trusted Python. They define a
`SCAFFOLD` and `GROUPS`; each group specifies its identity, fragment and states.
Inspect the distributed definitions in `config/starting_candidates.py` before
building. A read-only check is:

```bash
python -m redox build --config config/starting_candidates.py --only pmdi --dry-run
```

A successful dry run prints `[valid]` and creates no library. Normal construction
writes graft/conformer files and merges state rows into `library/manifest.csv`.
Calculation and analysis stages write separate `calcs/` and `results/` artifacts.
Do not equate a successful structure check with route or property qualification.

Downloaded datasets are described in `docs/DATASETS.md`; its dated inventories
are historical observations, not guarantees that those files exist in a clone.
New-source admission requires its own provenance and scientific checks.
