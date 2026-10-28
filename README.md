# Redox flow battery property pipeline

Computational screening of redox-active groups for post-functionalization of
preformed Merrifield-type chloromethyl-polystyrene, targeting **stability** and
**redox potential** in a nonaqueous (MeCN / PF₆⁻, Fc/Fc⁺-referenced) electrolyte.
Capped benzyl structures model individual graft sites; explicit finite oligomers
probe the local polymer environment. These model levels are reported separately.

This is the **refactored, self-contained** tree: maintained source, configuration,
curated input data, documentation, tests, and small cached examples. Bulk
calculation archives, campaign outputs, and manuscripts live in the original
`0_redox` workspace and are intentionally not carried over here.

## Layout

| Path | Contents |
| --- | --- |
| `src/redox/` | The installable package: build, QM, properties, screening, workflows, validation |
| `config/` | Candidate definitions, electrolyte/reference settings, execution policy |
| `data/raw/` | Curated inputs: candidate workbooks, two-wave MeCN benchmark |
| `docs/` | Modeling design, design axes, datasets, grafting notes |
| `scripts/` | Pipeline, mining, plotting and validation entry scripts |
| `tests/` | Pytest suite (self-contained via `pythonpath = ["src"]`) |
| `examples/` | Small cached thermochemistry example |
| `validation/nitro_prior/` | Source-bound Nitro force-repeatability evidence (see below) |

## Install

The package only links `redox`; the heavy environment (torch / fairchem /
gpu4pyscf / xtb) is built separately and lives in the `redox` conda env.

```bash
pip install -e .          # into the `redox` env, after the env is built
```

## Command-line entry points

`python -m redox --help` (or the installed `redox` command) lists every stage:

```
redox build        # build candidate structures (--dry-run validates only)
redox uma          # gas-phase warm starts (explicit backend qualification)
redox dft          # DFT/SMD with checked warm-start provenance
redox motifs       # plan or execute registered candidates
redox uma-check    # inspect UMA configuration without loading CUDA
redox status       # inspect all declared states and QC
redox review       # rebuild and audit tables in a fresh directory (no QM)
redox launch       # plan current DFT/UMA shards (--execute submits)
redox resources    # discover free GPUs through the system reservation gate
```

## UMA backend policy

The backend is chosen from the **actual reserved GPU's** compute capability, not
from a config preference:

- **Stock (legacy) backend by default** — runs on a V100 (cc 7.0) and any
  non-Nitro-capable card. It is deterministic (force repeatability ≈ 10⁻⁶ eV/Å).
- **Nitro** is selected only on a Nitro-capable card (CUDA compute capability
  8.x–12.x, i.e. an H100). It requires Hopper/Ampere-class hardware and **cannot
  run on a V100**.
- Nitro additionally stays blocked by a production-qualification gate
  (`config/uma_nitro_qualification.json`) until its force repeatability passes the
  1×10⁻⁴ eV/Å tolerance. The retained H100 evidence under
  `validation/nitro_prior/` currently shows it does **not** (fused-kernel TF32 +
  atomic scatter nondeterminism), so production falls back to the stock backend.

GPU work must go through the system reservation gate (`gpu_reserve`); the pipeline
verifies lease ownership before creating a CUDA context.

## Tests

```bash
pytest            # self-contained; uses this tree's src via pyproject pythonpath
```

## Documentation

- [Modeling the system](docs/MODELING.md)
- [Design axes and reliability](docs/DESIGN_AXES.md)
- [Datasets](docs/DATASETS.md)
- [Practical grafting](docs/PRACTICAL_GRAFTING.md)
- Agents should read [AGENTS.md](AGENTS.md) before changing or running the pipeline.
