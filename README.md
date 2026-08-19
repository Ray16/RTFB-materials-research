# Redox-active Merrifield monomers

Computational screening of redox-active functional groups grafted onto the Cl site of the
Merrifield-resin (chloromethyl-polystyrene) monomer, targeting **stability** and **redox
potential** in a nonaqueous electrolyte.

Operating instructions for Claude Code: [`CLAUDE.md`](./CLAUDE.md).

## System

- **Scaffold:** Merrifield monomer (chloromethyl-polystyrene); decorate the benzylic Cl
  site with redox-active groups (pyridine-derived / multi-electron systems, and others).
- **Electrolyte:** acetonitrile (ε ≈ 37.5), **SMD** implicit solvation, **PF₆⁻** counterion.
- **Redox potential:** `E° = −ΔG° / (nF)` (n = electrons transferred).
- **Reference dataset:** D3TaLES — https://d3tales.as.uky.edu/database/

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
| 2 | **DFT + SMD(MeCN)** geom-opt (PySCF) | **yes** | final solvated geometries → λ, RMSD, ΔG, E° |

- Warm-starting Tier 2 from Tier 1 cuts DFT optimization steps sharply.
- **Scaling tier:** if the library grows, insert **xtb + ALPB(MeCN)** solvated opt between
  Tiers 1 and 2, reserving DFT+SMD opt for a shortlist.
- Cross-check / open MLIP alternative: **AIMNet2** (organic, charge-aware) or `eSEN-omol`.
- Absolute E° from `E° = −ΔG/nF` is referenced to Fc/Fc⁺ (see `config/electrolyte.py`).

## Setup

`setup_env.sh` creates the `redox` env, installs pinned deps, and verifies. It
auto-detects GPU vs CPU; override with `CUDA_TAG`.

**GPU machine** (auto-detected → `cu124`):
```bash
./setup_env.sh
```
If your driver's CUDA version (`nvidia-smi`, top-right) is **below 12.4**, pick the
matching PyTorch wheel tag instead, e.g. `CUDA_TAG=cu118 ./setup_env.sh`.

**CPU-only machine** (auto-detected when no `nvidia-smi`, or force it):
```bash
CUDA_TAG=cpu ./setup_env.sh
```

Then activate and check:
```bash
source $(conda info --base)/etc/profile.d/conda.sh && conda activate redox
python scripts/check_env.py     # ENV OK; "cuda available: True" on GPU, False on CPU
```

UMA weights are gated — request access to `facebook/UMA` on HuggingFace, then
`huggingface-cli login`. Adding a pipeline dependency? Update `requirements.txt` +
`scripts/check_env.py` so the env stays reproducible.

## Layout

```
data/raw/          monomer structures, D3TaLES exports (inputs)
data/processed/    cleaned/derived datasets
config/            redox-group definitions, run params
src/redox/         build · uma · redox · descriptors
library/           generated decorated monomers (SDF/XYZ) + manifest
calcs/uma/         UMA optimization + energy outputs
calcs/dft/         DFT + SMD single-point outputs
results/           tables, plots, final descriptors
scripts/           runnable entry points
notebooks/         exploration
```

## Pipeline

1. `src/redox/build.py`  — decorate Cl site with redox groups → `library/`
2. `src/redox/uma.py`    — UMA relax neutral/ox/red states → `calcs/uma/`
3. `src/redox/redox.py`  — DFT+SMD single points → E°, ΔG, λ → `results/`
