# CLAUDE.md

Operating instructions for Claude Code in this repo. Project spec/background lives in
[`README.md`](./README.md) — read it for context, keep this file to instructions only.

## Environment
- Always activate the env first:
  `source ~/miniforge3/etc/profile.d/conda.sh && conda activate redox`
- Never install into `base`. Add packages to the `redox` env.
- Verify with `python scripts/check_env.py` after any env change.

## Guardrails
- Do **not** delete or overwrite anything in `data/`, `library/`, `calcs/`, or `results/`
  without explicit confirmation — calc outputs are expensive to regenerate.
- Do **not** commit model weights, large calc outputs, or HF tokens.
- Never hard-code or echo HuggingFace / API tokens into files or logs.
- Prefer open-source, reproducible tools (PySCF for DFT+SMD) over proprietary ones.

## Compute
- Pipeline: RDKit conformer → **UMA** (charge+spin) gas-phase pre-opt → **DFT+SMD(MeCN)**
  solvated geometry optimization for final structures/energies. Geometries are optimized
  **in solvent** (UMA is only the pre-optimizer; it cannot do solvation).
- Redox spans multiple charge states — every structure must carry an explicit **charge**
  and **spin multiplicity**; never assume neutral singlet.
- Solvent/referencing params live in `config/electrolyte.py` (acetonitrile, ε=37.5, SMD,
  PF₆⁻, Fc/Fc⁺ reference) — read from there, don't hard-code.
- One molecule = one stable ID; keep its calcs in a dir keyed by ID + redox state.
- GPUs: 4× V100 32GB — set `CUDA_VISIBLE_DEVICES` when launching parallel jobs.

## Conventions
- Generated structures → `library/`; raw inputs → `data/raw/`; don't mix.
- Config (redox groups, run params) → `config/`; reusable code → `src/redox/`.
- When adding a redox group, record its charge/spin per redox state in the config.

## Verification
- After code changes, run the relevant smoke test before reporting success.
- Report failures with the actual error; don't claim success on unverified steps.

---

## User instructions
<!-- Add project-specific instructions for Claude below. -->
