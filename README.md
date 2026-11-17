# Redox motifs for polymer grafting

Research software for assessing redox-active groups attached to preformed
Merrifield-type polymers. It constructs exact precursor/graft models, organizes
quantum-chemical calculations, and evaluates redox potentials, attachment
thermodynamics and finite oligomer responses with explicit method and quality
information.

This tree is the cleaned library. The manuscript, slides, dated campaign
directories and bulk calculation outputs remain in
`/nfs/lambda_stor_01/homes/rzhu/0_redox`. See [MIGRATION.md](MIGRATION.md) for
what moved, what stayed byte-identical, and what was left behind.

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

GPU work still goes through `gpu_reserve`. This refactor does not launch DFT
or Hessian jobs.

## Tests

```bash
OMP_NUM_THREADS=1 python -m pytest tests -q
```

The suite is the scientific-contract subset. Tests that only opened dated
campaign directories, the manuscript tree, or bulk `results/` were not copied.
Those records stay in the source repository.
