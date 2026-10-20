# Redox flow battery property pipeline

**Current results:** [results/README.md](results/README.md) indexes grafting free
energies, polymer models, candidate tables, validation figures and the slide deck.

Computational screening of redox-active functional groups grafted onto the Cl site of the
Merrifield-resin (chloromethyl-polystyrene) monomer, targeting **stability** and **redox
potential** in a nonaqueous electrolyte.

**Current correctness review:** rigid-motion contamination in native thermochemistry
has been corrected and independently checked. Current candidate exports use the
corrected thermal partition and matched Fc reference, with field-specific QC.
Earlier unprojected numerical tables are historical; broader scientific validation
remains incomplete. See the [correction evidence and status](results/validation/rigid_mode_thermal_2026-10-18/README.md).

Operating instructions for Claude Code: [`CLAUDE.md`](./CLAUDE.md).

## Start here

This repository contains reusable property calculations, candidate definitions, and
saved scientific evidence. Merrifield-grafted redox molecules are the current
screening application.

| Task | Start with |
| --- | --- |
| Navigate the code and outputs | [Repository map](docs/REPO_MAP.md) |
| Develop or debug the pipeline | [Development guide](docs/DEVELOPMENT.md) |
| Find current guidance and historical reports | [Documentation index](docs/README.md) |
| Understand each property and its QC | [Design axes](docs/DESIGN_AXES.md) |
| Review defects and refactoring decisions | [Code review](docs/CODE_REVIEW_2026-10-18.md) |
| Inspect the scientific continuation | [Correctness and expansion](docs/SCIENTIFIC_PROGRESS_2026-10-18.md) |
| Open the latest new-candidate results | [Candidate results](results/candidates/README.md) — workbook, CSV and matched comparisons |
| Inspect the completed full-property arithmetic review | [Reviewed property tables](results/validated/README.md) |
| Explore redox families and complementary datasets | [Candidate exploration](docs/CANDIDATE_EXPLORATION.md) |
| Choose an operational entry point | [Script guide](scripts/README.md) |

Agents should read [AGENTS.md](AGENTS.md) before changing or running the pipeline.
For checked tables, follow the snapshot in [LATEST_REVIEW.json](results/LATEST_REVIEW.json).
Root-level CSV exports may predate the current protocol. Individual module CLIs
still use those root paths by default; the isolated review workflow produces a
consistent export and records whether its inputs stayed unchanged.
For a first inspection, these commands read existing records:

```bash
export OMP_NUM_THREADS=1 OMP_THREAD_LIMIT=1 OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
PYTHONPATH=src python -m pytest -q
PYTHONPATH=src python -m redox --help
PYTHONPATH=src python -m redox status
PYTHONPATH=src python -m redox energy-audit
```

These caps keep local inspection within a small CPU budget. GPU launchers assign
their own disjoint CPU sets and native thread limits.

To rebuild tables in a fresh directory and run the independent arithmetic audit:

```bash
PYTHONPATH=src python -m redox review --out results/reviews/my_review
```

This writes an isolated review with stage logs and `review_status.json`. It performs
no quantum calculations and does not promote its tables to published results.
After an editable install, `redox --help` provides the same commands. Existing
`python -m redox.SUBPACKAGE.MODULE` entry points continue to work.

## System

- **Scaffold:** Merrifield monomer (chloromethyl-polystyrene); decorate the benzylic Cl
  site with redox-active groups (pyridine-derived / multi-electron systems, and others).
- **Electrolyte:** acetonitrile (ε ≈ 37.5), **SMD** implicit solvation, **PF₆⁻** counterion.
- **Redox potential:** `E° = −ΔG° / (nF)` (n = electrons transferred).
- **Reference data:** primary experimental MeCN/Fc measurements, D3TaLES and
  complementary computed collections; see the [source catalog](config/discovery_sources.json)
  and [experimental benchmark](data/raw/validation/two_wave_mecn/README.md).

### Descriptors
Reaction free energy · reorganization energy (λ) · structural RMSD between redox states ·
structural change on ox/red · electronic-property distributions · number of electrons
transferred · functional-group stability across redox states.

## Compute strategy (UMA-accelerated, solvated optimization)

Redox compares different **charge states**, so any MLIP used must accept charge (+ spin) —
this rules out charge-neutral universal potentials (MACE-MP, CHGNet, M3GNet, ORB, …).
Geometries are optimized **in the acetonitrile continuum** (SMD), since the solvated
structures set λ and the redox structural-change descriptors. UMA is gas-phase, so it
serves as the fast pre-optimizer that warm-starts the solvated DFT step — not the final
optimizer.

| Tier | Method | Solvent | Role |
|------|--------|:--:|------|
| 0 | RDKit ETKDGv3 ensemble → FF rank | no | conformer search → best seed |
| 1 | **UMA** (charge+spin, fairchem) | no | fast gas-phase pre-opt + gas-phase descriptors, all states |
| 2 | **DFT + SMD(MeCN)** — r2SCAN-D4/def2-SVP(D) opt // ωB97M-V/def2-TZVPD energy (one diffuse basis for every charge state and Fc), RI-J (PySCF / gpu4pyscf) | **yes** | final solvated geometries → λ, RMSD, ΔG, E° |

- Warm-starting Tier 2 from Tier 1 cuts DFT optimization steps sharply.
- **Scaling tier:** if the library grows, insert **xtb + ALPB(MeCN)** solvated opt between
  Tiers 1 and 2, reserving DFT+SMD opt for a shortlist.
- Cross-check / open MLIP alternative: **AIMNet2** (organic, charge-aware) or `eSEN-omol`.
- Absolute E° from `E° = −ΔG/nF` is referenced to Fc/Fc⁺ (see `config/electrolyte.py`).

## Setup

`setup_env.sh` creates the `redox` env, pins the validated QM libraries, and verifies. It
auto-detects GPU vs CPU; override with `CUDA_TAG`.

**GPU machine** (auto-detected → `cu128`, requires a compatible driver):
```bash
./setup_env.sh
```
Check the driver's supported CUDA version with `nvidia-smi` before installation.
`CUDA_TAG` must identify a compatible PyTorch 2.8 wheel. The recorded working
environment and native thermal version requirements are in the development guide.

**CPU-only machine** (auto-detected when no `nvidia-smi`, or force it):
```bash
CUDA_TAG=cpu ./setup_env.sh
```

Then activate and check:
```bash
source $(conda info --base)/etc/profile.d/conda.sh && conda activate redox
pip install -e . --no-deps       # editable install of the `redox` package (drops PYTHONPATH=src)
python scripts/check_env.py      # ENV OK; "cuda available: True" on GPU, False on CPU
```
`pip install -e .` is optional — every `python -m redox.*` command also works with a
`PYTHONPATH=src` prefix. Use `--no-deps` so pip never re-resolves the CUDA-matched torch stack.

UMA weights are gated — request access to `facebook/UMA` on HuggingFace, then
`huggingface-cli login`. Adding a pipeline dependency? Update `requirements.txt` +
`scripts/check_env.py` so the env stays reproducible.

## Layout

See [`docs/REPO_MAP.md`](./docs/REPO_MAP.md) for a full, up-to-date guide. In brief:

```
data/raw/          monomer structures, D3TaLES exports (inputs)
config/            redox-group defs, run params, electrolyte/referencing (+ project.json)
src/redox/         the `redox` package: core/ (paths, energy protocol), build/ (3D structures),
                   qm/ (UMA, DFT+SMD, single points), properties/ (E°, integrity, stability,
                   λ), screening/ (scorecard, Pareto), validation/ (benchmarks, D3TaLES)
library/           generated decorated monomers (XYZ) + manifest
calcs/uma/         UMA pre-optimization evidence
calcs/dft/         DFT + SMD geometry, addressed energies and thermal evidence
results/           property tables (CSV) + figures/{candidates,validation,reorg,pipeline}/
scripts/           entry points: pipeline/ (production chain), analysis/ (candidate
                   sensitivity), validation/{orop,solvation,reorg_d3tales,reorg_literature}/,
                   mining/ (D3TaLES), fleet/, polaris/, plotting/{candidates,validation,reorg,
                   pipeline}/ (mirrors results/figures/, one shared plot_style)
tests/             data-contract tests (pytest)
archive/           historical location; empty in the reviewed working tree
```

## Pipeline

New motifs can be added to `config/custom_candidates.py` and run with
`scripts/pipeline/run_motifs.py` (plan by default; `--execute --gpu INDEX` runs each compute
stage through the reservation gate). Input validation, a new-motif embedding test, a cached
integration run, and audited table/workbook exports are documented in
[`docs/VALIDATION_REVIEW_2026-10-10.md`](./docs/VALIDATION_REVIEW_2026-10-10.md).
That integration snapshot is historical. New candidate proposals live under
`results/discovery/`. Current repairs, calculations and evidence are identified by
[`results/LATEST_REVIEW.json`](results/LATEST_REVIEW.json) and
the [scientific progress guide](docs/SCIENTIFIC_PROGRESS_2026-10-18.md).

1. `src/redox/build/`      — decorate the site with redox groups → `library/`
2. `src/redox/qm/uma.py`    — UMA charge/spin pre-opt of every state → `calcs/uma/`
3. `src/redox/qm/dft.py`    — DFT+SMD geometry optimization in solvent → `calcs/dft/<id>/<state>/`
   (RI-J default; optional `--torsion-scan` for floppy species)
4. `src/redox/qm/sp.py`     — single points at the ONE active protocol (`redox.core.protocol`,
   ωB97M-V/def2-TZVPD for every state + Fc) → protocol-addressed `sp/<hash>.json` records
5. `src/redox/qm/native_thermal.py` — derivative-checked solution thermochemistry,
   with rigid motions excluded and explicit spin treatment; failed fields stay blank.
6. `src/redox/properties/` — `potentials` · `integrity` · `stability` · `reorg` · `lambda_outer`:
   E° vs live Fc, charged-state integrity, ΔG_disp, λ_O/λ_R (Nelsen), λ_o (PCM/Born) → `results/`
7. `src/redox/screening/` — `scorecard` → `pareto`: per-pool scorecard (contiguous path from the resting
   state, INCOMPLETE propagation) and the σ-aware Pareto shortlist

After calculations and QC, derived tables can be rebuilt in isolation with
`redox.workflows.tables`. `scripts/pipeline/finalize_after_dft.sh` is the explicit
publication chain: it consumes saved cross points and writes canonical results
and figures. It no longer starts quantum calculations. What each axis means and how accurate it is:
`docs/DESIGN_AXES.md`.

Reorganization energies are cross-checked against D3TaLES by reproducing **their exact level of
theory** — IP-tuned LC-ωHPBE/def2-SVP gas, using their per-molecule tuned ω — via
`scripts/validation/reorg_d3tales/validate_reorg_worker_d3tales.py` (results in
`results/d3tales_reorg_validation_d3level/`).
Matching their level reproduces their λ; our production level (ωB97M-V/def2-TZVPD, SMD, uniform diffuse)
differs by functional+basis, which is the more correct treatment for anions (see FINDINGS #10).
