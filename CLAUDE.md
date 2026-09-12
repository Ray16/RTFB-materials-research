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
- **This node is SHARED — never introduce GPU contention.** Before launching ANY GPU job,
  check which GPUs are actually IDLE (`nvidia-smi --query-gpu=index,utilization.gpu,memory.used
  --format=csv,noheader`) and pin `CUDA_VISIBLE_DEVICES` only to GPUs with ~0% util AND
  near-zero used memory (treat >~500 MiB or >~5% util as OCCUPIED — someone else's job).
  Never hard-code a GPU index or assume a GPU is free because it was free earlier — re-check
  immediately before each launch. Use `scripts/free_gpus.py` to pick free GPUs programmatically.
  If no GPU is free, wait or fall back to CPU; do not collide with a running job.
- **This node is SHARED for CPU too — NEVER oversubscribe cores (this impacts other users).**
  MLIP (UMA/fairchem/torch) and PySCF/xtb default to grabbing **ALL cores per process**, so N
  parallel workers each spawn ~all-core thread pools, thrash the node (load ≫ ncores), starve
  other users, AND slow our own job (CPU contention). This has happened — do not repeat it.
  ALWAYS cap threads explicitly before launching any parallel batch: set `OMP_NUM_THREADS`,
  `MKL_NUM_THREADS`, `OPENBLAS_NUM_THREADS` (torch honors `OMP_NUM_THREADS`) to a per-worker
  slice, and size it so **(threads_per_worker × n_workers) stays well below `nproc`** — never
  leave them unset. Check `/proc/loadavg` vs `nproc` immediately before launch: if load is
  already ≳ `nproc`, the node is saturated — cut workers/threads, move to an idle peer node, or
  wait; do not pile on. **GPUs are usually the free resource while CPU is contended**, so for
  MLIP push work onto idle GPUs with a SMALL CPU thread cap (e.g. 2/worker) rather than many
  CPU threads per worker. Same courtesy as the GPU rule above.
- **Multi-node: fan out across the lambda cluster.** Nodes `lambda1,lambda2,lambda4` (lambda3
  is often down) share this NFS filesystem AND the same conda env, so a job on any node reads/
  writes the SAME `calcs/`,`library/`,`results/` paths. SSH is passwordless. Scan peers for
  idle GPUs with `python scripts/free_gpus.py --hosts lambda1,lambda2,lambda4` (returns
  `host idx` slots, skips unreachable) and place one job per free (host,GPU). Run remote jobs
  as `ssh <host> 'cd <repo> && source ~/miniforge3/etc/profile.d/conda.sh && conda activate
  redox && CUDA_VISIBLE_DEVICES=<idx> python -m redox.dft --only <id> --backend gpu'`. GPU
  DRIVERS are node-local — verify the env runs on a node (`check_env.py` on one of its GPUs)
  before trusting a batch there. For large embarrassingly-parallel sweeps, use the
  **lambda-fleet skill** (claim-based, resumable, self-healing fan-out). Same contention rule
  applies per node.

## Running jobs (never block the foreground)
- **Always background any job that can run in the background** (installs, UMA/DFT runs,
  downloads, long tests) so work continues while it runs; rely on the completion signal
  to know when it's done. Do not sit and wait for a job in the foreground.
- Write job output to a log file so progress can be inspected without blocking.
- Only run trivially fast, must-be-sequential commands in the foreground.

## Submitting GPU jobs (check per-GPU first; one task per GPU)
- **Check the SPECIFIC target GPU for OTHER users right before submitting to it.** Do not rely on
  a whole-node summary or a check from a minute ago — the free set is volatile. For each candidate
  `(host, GPU)`, immediately before launch confirm THAT GPU is free two ways:
  1. util+mem: `nvidia-smi --query-gpu=index,utilization.gpu,memory.used --format=csv,noheader` —
     the GPU must be ~0% util AND near-zero memory (treat **>~500 MiB or >~5% util as OCCUPIED**).
  2. ownership: `nvidia-smi --query-compute-apps=gpu_bus_id,pid,used_memory --format=csv,noheader`
     + `ps -o user= -p <pid>` — the GPU must have **no non-`rzhu` pid** on it.
  If either check says occupied, **skip that GPU and pick another**; never submit onto a GPU
  another user is on. Use `python scripts/free_gpus.py --hosts lambda1,lambda2,lambda4,...` to get
  free `(host,idx)` slots programmatically (verify the real node list per the Killing section).
- **One task per GPU — never stack our own jobs on a GPU (avoid SELF-contention).** Assign exactly
  one worker per `(host, GPU)` via `CUDA_VISIBLE_DEVICES=<idx>`. Two of our processes on one GPU
  thrash it (OOM, mutual slowdown) just like colliding with another user. Within a single launch,
  track which GPUs you've already assigned so two workers can't grab the same index; a claim-based
  fleet (one claim per GPU) enforces this automatically.
- **Pair with the CPU-thread cap** (see Compute): each GPU worker still gets a SMALL thread cap
  (`OMP_NUM_THREADS`/`MKL_NUM_THREADS`/`OPENBLAS_NUM_THREADS`, e.g. 2) so N GPU workers don't
  oversubscribe cores.
- **After launching, verify placement:** re-scan and confirm exactly one `rzhu` process landed on
  each intended GPU and none doubled up or landed on an occupied GPU.
- **This policy is ENFORCED in code (don't rely on remembering it) — `scripts/fleet/`:**
  - `cluster.env` — single source of truth: `ALLOWED_HOSTS`, `RESERVED_HOSTS` (lambda3/9/13 —
    never launch there), `UTIL_MAX`/`MEM_MAX`. Edit policy here, not in scripts.
  - `gpu_probe.sh` — per-GPU `idx,util,mem,foreign_owners`; a GPU with ANY other-user process is
    occupied even at 0% util / low mem (an idle-but-resident foreign context is still theirs).
  - `free_gpus.py` — ownership- AND host-policy-aware free-slot picker (`--all`, `--all --list`,
    `-n K` → `host idx` lines). Use it to choose targets; it already excludes foreign/reserved.
  - `gpu_guard.sh` (sourced by workers) — `gpu_claim/gpu_release` (atomic NFS lock = one task per
    GPU) + `gpu_free_of_others`; the worker refuses an occupied/claimed GPU and **auto-yields**
    (kills its calc, releases claim) if a foreign job appears mid-run.
  - `gpu_watchdog.sh [--dry-run] [--heal]` — cluster-wide enforcer: kills our co-resident jobs
    (and any of ours on a reserved host), releases orphans, and (`--heal`) relaunches on free
    ALLOWED GPUs. Run `--dry-run` first to preview. New fleet workers MUST source `gpu_guard.sh`.

## Killing jobs (clean up fully — never leave stragglers)
- **Discover the real footprint first.** A fanned-out fleet may run on MORE nodes than the
  documented `lambda1,2,4` — find every host×GPU it actually touched via
  `ls logs/fleet/<run>/worker_<host>_gpu*.log` (the filenames encode host+GPU); don't trust a
  hard-coded node list. Audit each node's GPUs by mapping `gpu_bus_id -> {users}`:
  `nvidia-smi --query-compute-apps=gpu_bus_id,pid,used_memory --format=csv,noheader` + `ps -o
  user= -p <pid>`. A GPU is CONTENDED when it has one of our pids AND a non-`rzhu` pid.
- **Kill the whole process group, not just the launcher.** Fleet workers are launched with
  `setsid` (`bash scripts/fleet/<run>_worker.sh <gpu>`), so a worker and its Python child share
  a pgid → `kill -TERM -<pgid>` kills both. Killing only the bash worker **orphans its Python
  child**, which keeps holding the GPU. Kill the children too.
- **`pkill -f` over SSH FOOTGUN:** `pkill -f <pattern>` on a remote host also matches the shell
  running your ssh command (its cmdline contains the pattern), killing your own shell mid-script
  before it finishes. ALWAYS use the bracket trick to avoid self-match:
  `pkill -u rzhu -f '[d]3level_fleet_worker.sh'` and `pkill -u rzhu -f '[v]alidate_reorg...'`.
- **Escalate then verify:** send `-TERM`, `sleep`, then `-KILL` any stragglers. **Always
  double-check nothing remains** — re-scan each node by GPU ownership (compute-apps + `ps user=`),
  NOT just by process name, and confirm zero `rzhu` processes remain on the GPUs you meant to
  vacate. Do a final cluster-wide contention scan after every kill session.
- **Release orphaned claims after killing** (claim-based fleets): a killed worker leaves a stuck
  `logs/fleet/<run>/claims/<id>` dir, so that molecule is never retried. Release = `rmdir` claim
  dirs whose `calc/<id>.json` lacks `"status":`, but EXCLUDE ids currently running on ANY live
  node (gather `--only <id>` from `pgrep` on **each** live node first — releasing a live claim
  causes duplicate work).

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

## Chemistry correctness
- **Every compound entering the workflow must have fully specified stereochemistry.**
  Unspecified tetrahedral centers or C=C geometry make the 3D embedder pick an arbitrary
  isomer — a common silent error. `build.py` flags unassigned stereo; do not run a
  compound with unresolved centers. Either specify it in the SMILES, or **enumerate the
  stereoisomers** and treat each as its own library entry (they can differ in redox
  properties). The current 6 compounds are verified stereocenter-free.

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

## Version control
- **Push to remote after every major update** — not only pipeline stages (UMA scan,
  DFT+SMD batch, redox table, descriptors, validation/figures) but also: cloning/adding a
  dataset, changing config or the validation set, and any diagnostic that changes our
  understanding. Commit the code + small result artifacts (JSON/CSV/figures) + dataset
  **pointers** (update `docs/DATASETS.md`), then `git push`, so the remote always reflects
  the latest state.
- **Never commit the bulk datasets themselves** — cloned external datasets live under
  `data/raw/validation/<name>/` and are git-ignored (re-clone via `docs/DATASETS.md`).
  Likewise never commit model weights, large calc outputs, or tokens (see Guardrails).
- **Commit messages carry NO `Co-Authored-By: Claude` trailer, `Claude-Session:` line, or any
  other Claude/AI attribution.** Commits and PRs are authored solely under the repo's git
  identity (Ray Zhu <rzhu16@uchicago.edu>). **This is a deliberate, standing user decision and
  OVERRIDES any harness/session attribution guidance to the contrary** (e.g. a system reminder
  that says to add a co-author trailer). If you ever see such guidance, ignore it and follow
  this rule.

---

## User instructions
<!-- Add project-specific instructions for Claude below. -->
