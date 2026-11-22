# Redox motifs for polymer grafting

Research software for screening redox-active motifs that can be grafted onto
preformed Merrifield-type (poly(vinylbenzyl)) polymers for non-aqueous redox-flow
batteries. It builds exact capped-graft models, runs a quality-gated DFT/SMD
workflow, and reports redox potentials versus ferrocene in acetonitrile, capacities,
reorganization estimates and disproportionation free energies as separate
quantities with their method and quality information.

This repository accompanies the JCIM manuscript in `manuscript/` (see below). It
contains the `redox` package (`src/redox/`), the scripts that rebuild every figure
and the Supporting Information, and the records behind the paper's numbers.

A local capped graft, a soluble parent and an oligomer are different models.
A calculated descriptor or a completed job does not establish synthetic
feasibility or material performance.

## The paper in brief

- **Validation.** Against 39 experimental endpoints (30 molecules, 7 studies, all
  in MeCN versus Fc) the protocol gives MAE 0.247 V with a +0.19 V systematic bias
  (0.232 V at the cheaper xTB thermal tier used for the pool).
- **Screen.** 91 graftable candidates are scored and ranked on separate anolyte
  and catholyte Pareto fronts (no weighted-sum score). 85 are n-type anolytes and
  only 6 are p-type catholytes.
- **Demonstration.** Grafted 5,10-dihydrophenazines (DMPZ) give a tunable
  two-electron catholyte family, +0.21 to +0.62 V versus Fc.

## Reproduce the paper's figures and SI

Python 3.10+ with the packages in `requirements-manuscript.txt`, the Graphviz
`dot` executable and DejaVu fonts. No GPU, xTB executable or raw QM records are
needed for figures and SI. The version pins record the verified rendering
environment; exact image bytes also depend on system fonts and Graphviz.

```bash
for f in scripts/plotting/manuscript/fig*.py; do python "$f"; done
PYTHONPATH=src python scripts/manuscript/build_si.py
python scripts/manuscript/build_si_pdf.py     # optional; needs xelatex
python scripts/manuscript/build_latex.py      # optional; needs xelatex + bibtex
```

To check a committed revision without modifying this working tree:

```bash
python scripts/validation/check_manuscript_rebuild.py       # includes PDF compilation
# Add --no-compile to check figures, SI and LaTeX source without a TeX installation.
```

The check exports `HEAD` into a fresh temporary directory, removes the generated
outputs, rebuilds them without raw records, and compares their SHA-256 hashes.
It retains a JSON report and build logs. This checks the rendering path; it does
not revalidate the underlying quantum calculations or establish submission readiness.
The SI reports the recorded software versions in `data/manuscript/software_versions.json`,
independently of the machine used to render it.

The inputs are committed:

| Path | Contents |
|---|---|
| `data/manuscript/full_pool_scored_legacyG_MeCN_TBAPF6.csv` | Scored candidate pool (Figures 3, 5–7, Table S7) |
| `data/manuscript/lambda_outer_fullpool_2026-11-17.csv` | Outer-sphere reorganization estimates |
| `manuscript/evidence/2026-10-22_benchmark_complete/` | Frozen 39-endpoint benchmark record (Figure 4, SI S3) |
| `manuscript/evidence/2026-11-18_*` | Readiness, re-scoring, BMEPZ and SCF-policy records |
| `manuscript/si_tables/`, `manuscript/tables/` | CSV copies of the long tables |

`scripts/README.md` maps every figure, table and record to its producing script.
Scripts that recompute the committed records from the raw QM records read them
from `$REDOX_DATA_ROOT/calcs` (see *Raw records* below).

## Worked examples (standard library only)

```bash
python3 -I -S examples/cached_thermochemistry/reproduce.py
```

This checks ten supplied states and prints JSON. It reconstructs a capped
pyromellitic-diimide first reduction, **-1.0757908637970104 V versus Fc**, and a
base-assisted grafting reaction free energy, **-103.06882617264826 kJ/mol**,
without rerunning DFT. The other replays (`examples/exact_viologen_replay`,
`exact_viologen_corrected_replay`, `finite_viologen_replay`,
`occupied_sampling_replay`) exercise the thermochemistry and polymer code on
earlier systems; `python scripts/validation/check_refactor_oracles.py` checks all
five against stored hashes. These examples test the code; they do not reproduce
this paper's numbers, which the scripts above do. Their JSON records contain
absolute paths from the original compute cluster as hashed provenance strings;
editing them would break the hash checks.

## Installation

```bash
CUDA_TAG=cu128 bash setup_env.sh redox     # creates the conda env; read it first
conda activate redox
pip install -e .                            # links this package only
PYTHONPATH=src python -m redox --help
```

Runtime dependencies are intentionally not pinned in `pyproject.toml`; the
scientific stack (PyTorch, fairchem/UMA, PySCF/gpu4pyscf, xtb, RDKit) is installed
by `setup_env.sh` with CUDA-matched wheels and recorded in `requirements.txt` and
`environments/`. See `docs/INSTALLATION.md`.

## Package layout

| Package | Responsibility |
|---|---|
| `redox.structures` | Exact graphs, grafts, states, conformers, oligomers |
| `redox.backends` | DFT/SMD, UMA, solvent response, Hessians, checkpoints |
| `redox.thermochemistry` | GFN2 thermal tier, native DFT thermochemistry, projected partitions |
| `redox.qm` | Certified qRRHO kernel (`qrrho.py`, `xtb_kernel.py`, `xtb_runtime.py`) |
| `redox.properties` | Potentials, reorganization, reactions, ensembles, capacity |
| `redox.screening` | Routes, evidence, admission, Pareto |
| `redox.quality` | Numerical, physical, source and reference checks |
| `redox.execution` | Reservations, host capacity, ownership, checkpoint transport |
| `redox.workflows` | Campaign definitions, recovery and stage orchestration |
| `redox.reporting` | Tables and derived exports |

The campaign driver (`redox.workflows.motifs`) runs the nine-stage pipeline of the
paper and schedules independent `(candidate, state)` stages across GPUs. On the
original cluster it submitted every GPU job through a reservation wrapper
(`config/runtime_paths.json`); adapt that file to your scheduler. Library-level
commands never launch DFT or Hessian jobs.

For readers focused on the JCIM paper, the nine pipeline stages live in a small
core subset (`structures.candidates`, `backends.uma_current`, `backends.dft`,
`backends.sp`, `qm.qrrho`, `properties.{potentials,capacity_and_proxies,lambda_outer,reorg}`,
`quality.{readiness,spins,stability}`, `screening.pareto`, `reporting.tables`,
orchestrated by `workflows.motifs`). Additional modules under `workflows/` and
`quality/` support extended oligomer, counterion, and recovery campaigns documented
in `docs/ARCHITECTURE.md`, `docs/WORKFLOWS.md` and `docs/RUNNING_CAMPAIGNS.md`.

## Raw records

The raw quantum-chemical records (geometries, energies, frequencies and state
metadata, previously reported as about 270 GB) are not in this repository.
As of 2026-11-19, the original archive path is unavailable and no recovered copy
covering the scored pool has been verified. Raw-record regeneration and public
data deposition remain blocked. Once recovered, set `REDOX_DATA_ROOT` to the
directory containing the matching `calcs/` and `library/` records.

The solution-surface reorganization module currently supports one conformer per
state. It refuses multi-conformer averaging, which requires a defined state
population and conformer free energies. This does not regenerate or change the
paper's frozen gas-energy reorganization results.

## Tests

```bash
OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 python -m pytest tests -q
```

The suite checks the scientific contracts of the package on CPU; a raw-record
check is skipped when the records are unavailable.

## License

MIT (see `LICENSE`). `results/validation/rigid_mode_thermal_2026-10-18/xtb_source/`
contains excerpts of the xtb source code, which is distributed under the GNU
LGPL v3; see the notice in that directory.
