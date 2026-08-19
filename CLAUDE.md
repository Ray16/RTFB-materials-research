# CLAUDE.md

Operating instructions for Claude Code in this repo. Project spec/background lives in
[`README.md`](./README.md) — read it for context, keep this file to instructions only.

## Environment
- Always activate the env first:
  `source ~/miniforge3/etc/profile.d/conda.sh && conda activate redox`
- Never install into `base`. Add packages to the `redox` env.
- **Keep the env reproducible:** whenever you add a new Python package to the main
  pipeline, immediately add it to `requirements.txt` (and `setup_env.sh` if it needs
  special handling, e.g. a CUDA-matched index) so `./setup_env.sh` reproduces the env
  after a fresh `git pull`. Also add its import to `scripts/check_env.py`.
- torch must match the GPU driver's CUDA version — install it via `setup_env.sh` (pinned
  before other deps), not ad hoc, or fairchem's resolver will pull a mismatched build.
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
- Do not assume specific hardware — **detect** available GPUs/CPUs at runtime (e.g.
  `torch.cuda.device_count()`, `os.cpu_count()`) and scale to what's present.

## Running jobs (never block the foreground)
- **Always background any job that can run in the background** (installs, UMA/DFT runs,
  downloads, long tests) so work continues while it runs; rely on the completion signal
  to know when it's done. Do not sit and wait for a job in the foreground.
- Write job output to a log file so progress can be inspected without blocking.
- Only run trivially fast, must-be-sequential commands in the foreground.

## Parallelization (always maximize)
- The redox **states/molecules are independent** — treat every stage as embarrassingly
  parallel; never run them serially when they can fan out.
- **GPUs:** detect how many are present and distribute independent jobs one-per-GPU via
  `CUDA_VISIBLE_DEVICES`, keeping all of them busy. For UMA, batch structures through the
  calculator rather than looping. Fall back to CPU when no GPU is present.
- **CPUs:** use RDKit `numThreads=0` (all cores) for conformer generation; set
  `OMP_NUM_THREADS`/`MKL_NUM_THREADS` for PySCF/xtb and shard jobs across cores.
- **Avoid oversubscription:** partition cores/GPUs across workers; don't let N jobs each
  grab all threads. Prefer a job queue that pins each task to a GPU + a core slice.
- Make batch scripts resumable (skip states whose output already exists) so parallel
  reruns don't redo finished work.

## Conventions
- Generated structures → `library/`; raw inputs → `data/raw/`; don't mix.
- Config (redox groups, run params) → `config/`; reusable code → `src/redox/`.
- When adding a redox group, record its charge/spin per redox state in the config.

## Task tracking
- Keep `TODO.md` current: mark tasks done as you finish them, and add newly discovered
  tasks. It's the living worklist; `docs/PLAN.md` is the strategic phased plan.

## Verification
- After code changes, run the relevant smoke test before reporting success.
- Report failures with the actual error; don't claim success on unverified steps.

---

## User instructions
<!-- Add project-specific instructions for Claude below. -->
