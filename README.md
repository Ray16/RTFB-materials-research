# Redox motifs for polymer grafting

Research software for assessing redox-active groups attached to preformed
Merrifield-type polymers. It constructs exact precursor/graft models, organizes
quantum-chemical calculations, and evaluates redox potentials, attachment
thermodynamics and finite oligomer responses with explicit method and quality
information.

A local capped graft, a soluble parent and an oligomer are different models.
The software keeps their structures, charge/spin states, solvents, reference
conventions and validation status explicit. A calculated descriptor or a
completed job alone does not establish synthetic feasibility or material performance.

## Quick start: reproduce a cached result

From the repository root, Python 3.10 or later is sufficient:

```bash
python3 -I -S examples/cached_thermochemistry/reproduce.py
```

This standard-library example checks ten supplied states and prints JSON. It
reconstructs the capped PMDI-ethyl first reduction, approximately **−1.075791 V
versus Fc**, and the defined DBU separated-ion reaction, approximately
**−103.068826 kJ/mol**. It needs no GPU, model download or original data paths.
Use `--out result.json` to save the JSON. It checks cached arithmetic and
provenance; it does not rerun DFT or establish experimental accuracy.
[Example inputs and scientific scope](examples/cached_thermochemistry/README.md).

## Install and inspect

The package requires Python ≥3.10. Python 3.11 is used in the working scientific
environment. Use an editable checkout: configuration and scientific records live
beside `src/`, rather than inside a standalone wheel.

```bash
python -m pip install --no-deps --no-build-isolation -e .
python -m redox --help
python -m redox workflows
python -m redox workflows --category controls --json
```

Editable installation requires setuptools and pip. These inspection commands
use the standard library. The package deliberately declares no automatic runtime
dependencies, so installing it does **not** install the scientific environment.
See [installation and dependencies](docs/INSTALLATION.md) before building
structures or running calculations.

## Workflow

| Stage | Inputs and operation | Outputs and admission |
| --- | --- | --- |
| Design and route review | Exact precursor/product graphs, attachment bond and source evidence | Candidate definitions and proceed/hold/exclude route decisions |
| Structure construction | Trusted Python configuration with scaffold, fragments and states | `library/manifest.csv`, graft identities and conformer XYZ files |
| Geometry and electronic calculations | Declared UMA/DFT methods, solvent, charge/spin and resource policy | Geometry/checkpoint files and protocol-addressed energy records in `calcs/` |
| Thermochemistry and property checks | Each state's own modes/thermal terms and compatible references | Qualified redox/reaction/oligomer quantities; missing components remain pending |
| Publication and replay | Frozen evidence and explicit property definitions | Reviewed tables, figures, manuscript and portable examples |

Inspect a command's `--help` before using it. `redox motifs` plans a registered
campaign by default; `--execute` authorizes execution. Production orchestration,
calculation stages, numerical controls and recovery tools have different roles:
[maintained workflow guide](docs/WORKFLOWS.md). Detailed commands, previous
calibration notes and measured campaign timings are preserved in
[running campaigns](docs/RUNNING_CAMPAIGNS.md).

Full molecular DFT/Hessian jobs require the validated 32-GB GPU class and the
system `gpu_reserve` gate. Native thread pools, CPU allocation, host RAM and
source-bound QC remain mandatory. UMA backend eligibility depends on the actual
reserved GPU and its qualification evidence. Installation on another host does
not establish numerical reproducibility there.

## Results, manuscript and reproducibility

- [Results index](results/README.md): current outputs, field-specific QC, raw paths
  and source manifests; these are project runtime artifacts, not guaranteed inputs
  of a fresh clone.
- [Live manuscript](manuscript/overleaf/README.md): the shared Overleaf document.
  Frozen source exports and earlier manuscript versions retain their own scope.
- [Publication readiness](docs/PUBLICATION_READINESS.md): remaining scientific and
  release requirements. A permanent deposition identifier is still pending.
- [Modeling assumptions](docs/MODELING.md), [property definitions](docs/DESIGN_AXES.md)
  and [practical grafting](docs/PRACTICAL_GRAFTING.md).

The 39-endpoint parent benchmark is retrospective. Finite oligomer comparisons
retain their sampled configurations and ion conventions. Occupied ensemble
convergence, salt speciation, experimental exact-graft accuracy and material
accessibility have separate requirements. Frozen numerical replay and fresh
quantum-chemical reproduction are distinct validation scopes.

## Repository map

| Path | Purpose |
| --- | --- |
| `src/redox/` | Reusable structure, QM, property, workflow and validation modules |
| `config/` | Chemical definitions, physical protocols and execution policies |
| `data/` | Curated inputs and declared reference/source records |
| `examples/` | Portable cached numerical replays with scoped environments |
| `tests/` | Regression, scientific-contract and input-validation tests |
| `scripts/` | Campaign, validation, release and plotting entry points |
| `docs/` | Methods, installation, workflows and release limitations |
| `manuscript/` and `slides/` | Paper collaboration, evidence and presentation files |
| `library/`, `calcs/`, `results/` | Generated scientific records and local runtime outputs |

All new work belongs under this canonical project directory. Historical sibling
paths temporarily resolve through the relocation map in
`validation/workspace_relocation.json`; the running polymer campaign moves only
after its calculations finish. Do not rewrite immutable evidence paths to make
an incomplete relocation appear complete.

## Development and citation

See [contributing](CONTRIBUTING.md) for the lightweight checks and the procedure
for changing scientific behavior. Full tests use the declared scientific
runtime; a small standard-library/pytest subset is suitable for CPU CI.

The repository has an existing [MIT license](LICENSE). Third-party data, papers
and model weights retain their own terms. The manuscript is a working draft;
there is no final publication DOI to cite yet. For now, identify the repository
URL, exact commit and the specific calculation/evidence snapshot used. No final
publication metadata is implied by this README.
