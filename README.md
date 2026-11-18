# Redox motifs for polymer grafting

Research software for assessing redox-active groups attached to preformed
Merrifield-type polymers. It constructs exact precursor/graft models, organizes
quantum-chemical calculations, and evaluates redox potentials, attachment
thermodynamics and finite oligomer responses with explicit method and quality
information.

This tree is the single home for active work: the library (`src/redox/`),
campaigns (`results/`), and the manuscript (`manuscript/`). Historical slides,
dated campaign directories and bulk calculation outputs live in the read-only
archive at `/nfs/lambda_stor_01/homes/rzhu/archive/0_redox`. See
[MIGRATION.md](MIGRATION.md) for what moved, what stayed byte-identical, and what
was left behind.

A local capped graft, a soluble parent and an oligomer are different models.
A calculated descriptor or a completed job does not establish synthetic
feasibility or material performance.

## Reproduce a cached result

Python 3.10 or later, standard library only:

```bash
python3 -I -S examples/cached_thermochemistry/reproduce.py
```

This checks ten supplied states and prints JSON. It reconstructs the capped
PMDI-ethyl first reduction, **-1.0757908637970104 V versus Fc**, and the
defined DBU separated-ion reaction, **-103.06882617264826 kJ/mol**. It does
not rerun DFT.

The other portable replays are `examples/exact_viologen_replay`,
`examples/exact_viologen_corrected_replay`, `examples/finite_viologen_replay`
and `examples/occupied_sampling_replay`. The occupied-sampling replay imports
this checkout.

## Package layout

| Package | Responsibility |
|---|---|
| `redox.structures` | Exact graphs, grafts, states, conformers, oligomers |
| `redox.backends` | DFT/SMD, UMA, solvent response, Hessians, checkpoints |
| `redox.thermochemistry` | GFN2 diagnostic, native DFT thermochemistry, projected partitions |
| `redox.qm` | Certified qRRHO kernel (`qrrho.py`, `xtb_kernel.py`, `xtb_runtime.py`) |
| `redox.properties` | Potentials, reorganization, reactions, ensembles, capacity |
| `redox.screening` | Routes, evidence, admission, Pareto |
| `redox.quality` | Numerical, physical, source and reference checks |
| `redox.execution` | Reservations, host capacity, ownership, checkpoint transport |
| `redox.workflows` | Campaign definitions, recovery and stage orchestration |
| `redox.reporting` | Tables and derived exports |

`redox` command names are unchanged. `pip install -e .` links only this
package; runtime dependencies stay unpinned. Install into a scientific
environment built by `setup_env.sh`, and do not point that install at a
different checkout by accident.

```bash
PYTHONPATH=src python -m redox --help
```

GPU work still goes through `gpu_reserve`. Library-level commands do not launch
DFT or Hessian jobs; the campaign driver (`redox.workflows.motifs`) submits them
through the reservation gate and now schedules independent `(candidate, state)`
stages as a bag-of-tasks DAG across the reserved GPUs.

## Manuscript

The manuscript (method-led JCIM paper) is in `manuscript/`: the text
`manuscript/manuscript_draft_v0.md`, the journal (achemso) version in
`manuscript/jcim/` (built by `scripts/manuscript/build_latex.py`), the Supporting
Information `manuscript/supporting_information.md` (built by
`scripts/manuscript/build_si.py`), and the records behind specific statements in
`manuscript/evidence/`. The validation anchor is **MAE 0.247 V, +0.19 V bias over
39 endpoints** (30 molecules, 7 studies; 0.232 V at the pool's xTB thermal tier). Engineering plans behind recent pipeline changes are in
`docs/EXECUTION_REFACTOR_PLAN.md` and `docs/EFFICIENCY_SURVEY.md`.

## Tests

```bash
OMP_NUM_THREADS=1 python -m pytest tests -q
```

The suite is the scientific-contract subset. Tests that only opened dated
campaign directories, the manuscript tree, or bulk `results/` were not copied.
Those records stay in the source repository.

For a byte-exact comparison of all five cached replays with the original source:

```bash
python scripts/validation/check_refactor_oracles.py
```

Use the scientific Python environment for NumPy/RDKit-dependent replays. See
`baseline/refinement_oracles.jsonl` and `baseline/refinement_tests.txt` for the
latest refinement checks and `MIGRATION.md` for remaining architecture work.
