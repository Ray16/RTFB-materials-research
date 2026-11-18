# Plan: parallelize the per-candidate DFT screening campaign (`motifs.py`)

Repo: `/nfs/lambda_stor_01/homes/rzhu/0_redox_refactor`
Scope: scheduling + reporting only. No change to any QM/thermo kernel, protocol, or certificate.
Status: WRITTEN PLAN — no code edits, no GPU/DFT runs.

---

## 0. What exists today (confirmed by reading the source)

- `src/redox/workflows/motifs.py`
  - `commands()` (L118-141) builds a FLAT per-candidate list:
    `build → uma → dft → sp → [native_thermal → fc_reference_sync] → coverage → cross_points` per id,
    then ONE `review_tables` (L139) for the whole pool.
  - `main()` loop (L296 `for index,(label,compute,argv) in enumerate(jobs)`) runs them STRICTLY SERIALLY.
  - Only `dft` fans out: when `workers>1` (L299) it calls `run_parallel(state_jobs(...), zip(gpus,cpu_blocks))`
    (L89-116). `run_parallel` does `list(pool.map(one, jobs))` — a **barrier**: it waits for all 3 states of
    the *current* candidate before returning, so the slowest state (ox1) gates, and `sp`/next candidate wait.
  - GPU wrapping is already correct: `wrap()` (L291-293) prepends `taskset --cpu-list … gpu_reserve run <gpu> -- …`.
    `run_stage()` (L65-86) retries only gate refusals. Per-gid exclusive `flock` in `logs/pipeline_locks/<gid>.lock`
    (L259-271). Manifest written incrementally via `atomic_write_json(mp, manifest)` after every step (L304/315/318).
- `src/redox/workflows/stages.py` — `execute_graph(stages, run, *, workers, started, finished)` (L23-62):
  dependency-aware DAG scheduler with `workers` slots. Calls `run(stage, slot)` once per node after all its
  `stage.dependencies` succeed; `slot ∈ [0,workers)` is the disjoint-resource index; stops scheduling new work
  after the first failure but lets in-flight nodes finish. `validate_graph` (L9-20) rejects cycles/missing deps.
  **motifs.py does not use it; `reporting/tables.py` already does.**
- `src/redox/reporting/tables.py` — `review_tables(dest, workers=N)` (L129-228) ALREADY drives `STAGES`
  (a `TableStage` DAG, L40-68) through `execute_graph` when `workers>1` (L207-211), each stage in a fresh
  interpreter via `run(stage,slot)` (L177-185) with `taskset` affinity. `snapshot_inputs()` (L69-104) produces a
  per-file sha256 inventory saved as `input_inventory.json`; `review_tables` already diffs initial vs final inputs
  (L217-220 `changed_inputs`). It is a standalone module: `python -m redox.reporting.tables --out …`.
- Fleet/claim machinery:
  - `execution/chunked.py` — atomic NFS pending→worker `claim()` (L284-299), resumable `run_worker` (L418-488),
    `recover_local` (L491-524). **Explicitly CPU-only** (`module_for_stage` in `cpu_chunks.py` L21-31 raises
    "quantum jobs require their GPU workflow").
  - `execution/prepared_namespace_queue.py` — long-lived controller (`run`, L37-113) that observes/dispatches
    **already-prepared immutable namespaces** across the fleet via `charged_native_followup.dispatch_rows`;
    used by `route_molecular_panel.queue()` (L109-116). Needs the "prepared namespace" abstraction + per-row
    `namespace_sha256`; it is heavyweight and assumes a frozen roster.
  - `execution/launch.py` + `execution/gpu.py` — `launch.py` fans **static `--shard N:i`** DFT/UMA shards across
    all free GPUs of ONE host; each shard runs through `execution/gpu.py` which does the full
    `gpu_reserve run + preflight + CPU affinity + env` for one GPU. Backends `dft`/`uma` already accept
    `--shard`/`--all` (L28-42 of launch.py).

### Measured per-candidate (from the two `run_manifest.json` files)
build 0.3 s · uma 27-290 s (1 GPU) · dft 26-40 min on 3 GPUs but **gated by ox1**; neu/ox2 finish 10-26 min →
2 GPUs idle ~15 min · sp ~15 s (1 GPU) · reorg 5-6 min (1 GPU) · coverage 0.2 s · review_tables **46 min, serial,
and it FAILED (rc=1)** in the seed run after all QM succeeded.

---

## (a) Serial loop → bag-of-tasks GPU queue

### Decision: reuse `execute_graph` (stages.py) as the campaign scheduler. Do NOT reinvent, do NOT use the
### chunked/prepared-namespace queues for the single-node case.

Rationale:
- `execute_graph` already gives exactly the two things we need: (i) **intra-candidate ordering** via per-node
  `dependencies`, and (ii) **inter-candidate pipelining** via `workers` slots that start the next ready node the
  instant a slot frees — so a GPU freed by candidate A's `ox2` immediately picks up candidate B's `neu`/`uma`/`sp`,
  never idling on ox1 or on table-building.
- The GPU-safe wrapping (`wrap`/`run_stage`/`gpu_reserve`) and per-gid locks already exist in motifs.py and are
  reused verbatim inside the `run` callback.
- `chunked.py` is CPU-only by contract; `prepared_namespace_queue` requires pre-frozen namespaces and is aimed at
  the multi-node fleet — overkill and a poor fit for a mixed DAG (build/uma/dft/sp/native/coverage/reorg) that
  mutates `calcs/` in place. Keep it in reserve for the >1-node case (see "6 GPUs" below).

### New glue (all inside `motifs.py`; keep `commands()` as the source of the exact argv per label)

1. **`build_graph(ids, …) -> list[MotifNode]`** — replaces the flat `jobs` list. One node per unit of work:
   - `build:<id>` (CPU, no dep)
   - `uma:<id>` (GPU, dep: `build:<id>`)
   - `dft:<id>:<state>` for each state in `declared_rows()` for that id (GPU, dep: `uma:<id>`) — **this is the
     key change: per-(candidate,state) nodes replace the single `dft` barrier node.** Reuse `state_jobs()`
     (L89-95) to produce each state's argv (`--only <gid>:<state>`), which already targets the independent
     `calcs/dft/<id>/<state>/` directory.
   - `sp:<id>` (GPU, dep: all three `dft:<id>:*`)
   - `native_thermal:<id>` (GPU, dep: `sp:<id>`) — native mode only
   - `fc_sync` — **single shared node**, dep: every `native_thermal:*` (see correctness §c; must NOT be per-id concurrent)
   - `coverage:<id>` (CPU, dep: `sp:<id>`; in native mode also `fc_sync`)
   - `reorg:<id>` (`cross_points`, GPU, dep: `coverage:<id>` + `sp:<id>`)
   Node = small dataclass carrying `name`, `dependencies` (tuple), `label`, `argv`, `compute` (GPU bool).
   `validate_graph` is run first (free cycle/dep check, already in stages.py).

2. **`run(node, slot)` callback** — maps `slot → (gpus[slot], cpu_blocks[slot])` (the static per-slot GPU+CPU
   binding already computed at L243), then:
   - `argv2 = wrap(node.argv, gpus[slot], cpu_blocks[slot], node.compute)`
   - `proc, attempts = run_stage(argv2, log, env)` (verbatim reuse — keeps gate-refusal retry + provenance).
   - returns `proc.returncode` (the int code `execute_graph` expects).
   `workers = len(gpus)`, so one slot == one reserved GPU == one task per GPU (satisfies the one-task-per-GPU rule;
   `gpu_reserve run <gpus[slot]>` still arbitrates). CPU-only nodes (build/coverage/fc_sync) also take a slot but
   skip the gate (`compute=False` path in `wrap`) — acceptable, they are sub-second except fc_sync.

3. **Provenance**: keep the manifest exactly as today but populate it from `execute_graph`'s `started`/`finished`
   callbacks (same pattern tables.py uses, L169-195). On `started(node)` append
   `dict(stage=node.label, node=node.name, argv=…, log=…, started=…)` to `manifest["steps"]` and
   `atomic_write_json(mp, manifest)`; on `finished(node, code)` record `returncode/finished/attempts`. This keeps
   the `steps[].argv` gpu_reserve pattern and per-step logs byte-compatible with existing manifests/readers.

4. **Resumability**: today resume re-runs `commands()` into a fresh `resume_<ts>/` dir and re-executes every step
   (steps that find valid outputs are cheap/no-ops). Preserve that, and additionally add a per-node **`done(node)`
   predicate** so finished units are skipped without a subprocess: for `dft:<id>:<state>` check the existing
   state readiness used elsewhere (`redox.quality.readiness.inspect_states` / `read_result(id,state)` already
   imported at L320/329); for `reorg`/`sp` check their output artifacts. In `build_graph`, seed `execute_graph`'s
   `completed` set by pre-filtering nodes whose `done()` is true (or short-circuit inside `run`). This makes
   resume O(remaining work) instead of O(all candidates) and is safe because each state's `calcs/dft/<id>/<state>/`
   is already isolated. The `--resume` input/method guard (L200-214) is unchanged.

5. **CLI**: `--gpu 0,1,2` semantics unchanged; `workers=len(gpus)` as today (L218). The CPU preflight at
   L217-243 (host_capacity, cpu_allocation, per-worker `cpu_blocks`) is reused unchanged — it already sizes for
   `threads*workers` disjoint blocks, which is exactly what the DAG slots need.

### What is explicitly NOT rebuilt
- No new claim/journal layer, no SQLite, no NFS rename protocol for the single-node path — `execute_graph` +
  the existing per-gid `flock` + `gpu_reserve` already give safe, resumable, owner-aware dispatch on one node.

---

## (b) Decouple + incrementalize `review_tables`

Two independent moves:

### b1. Take it OFF the campaign critical path (low risk, do first)
`motifs.py` already has the mechanism: `--skip-review` (L151, L192-193) drops `review_tables` **and**
`fc_reference_sync`, and the `a.skip_review` branch (L332-348) instead does a cheap per-record verification
(`fc_live()`, `read_result`, `free_energy`, `inspect_states`) and ends the campaign with
`status="selected_outputs_verified"`. Plan:
- Make the campaign **default to the skip-review completion path** (verify selected records, write manifest,
  stop) and run `review_tables` as a **separate downstream job**: `python -m redox.reporting.tables --out
  results/<run>/tables --workers N`. (Keep `fc_sync` in the DAG per §c — only the *table build* is removed from
  the campaign; the Fc reference is still set.) This immediately removes the **46-min serial SPOF tail**: the
  seed-run `review_tables` failure (rc=1) would no longer discard a campaign whose QM all succeeded.

### b2. Parallelize the table build (already supported)
`review_tables(workers=N)` already runs the `STAGES` DAG through `execute_graph` (L207-211). Just invoke with
`--workers 8`. The DAG (L40-68) has wide independent branches (potentials, orop_benchmark, reorganization,
integrity, stability, capacity, outer_sphere all depend on nothing; only scorecard/pareto/workbook are the
narrow tail), so 46 min → roughly the longest single stage + the serial scorecard→pareto→workbook tail
(expect ~8-12 min at workers≈8). This is off the critical path, so even its wall time no longer blocks GPUs.

### b3. Incrementalize (only recompute changed candidates) — larger, staged work
`snapshot_inputs()` (L69-104) already fingerprints every per-record input (sha256) and `review_tables` already
persists it as `input_inventory.json` and computes `changed_inputs` vs the final snapshot (L217-220). Reuse that
as the incremental key:
- Add a persistent `results/tables_cache/` keyed by `input_inventory.json` entries. On a new review, diff the
  current `snapshot_inputs(source)` against the **previous review's `input_inventory.json`**; map changed files →
  affected candidate ids (paths under `calcs/.../<id>/...`).
- For the **per-record stages** (potentials, reorganization, integrity, stability, orop_benchmark — the ones whose
  cost scales with pool size), recompute only the changed ids and reuse cached per-record rows for the rest,
  keyed on the record sha256. This needs each of those stage modules to accept a "subset + cache" path; that is
  new per-stage glue, so stage it after b1/b2.
- The **pool-level aggregators** (scorecard, pareto, workbook) inherently need all rows but are cheap; they re-run
  over `cached ∪ recomputed`. Because tables.py already runs every stage in a fresh interpreter and fingerprints
  all inputs, the sha256-keyed cache preserves the independent-audit guarantee (a changed record → changed hash →
  forced recompute).
- Lowest-risk first increment: cache + skip at the **candidate granularity for the dominant stage(s) only**
  (profile first — likely `reorganization`/`potentials`), leaving everything else full-rebuild.

---

## (c) Correctness / ordering constraints that MUST be preserved

1. **Intra-candidate DAG edges** (encode exactly, do not relax):
   `uma → dft[neu,ox1,ox2] → sp → native_thermal`; `sp → coverage → reorg` (coverage GATES reorg);
   `reorg` also needs `sp`. These match the current flat order (dft before sp, coverage before cross_points,
   native_thermal after sp).
2. **`fc_reference_sync` is GLOBAL shared state and the single biggest hazard.** `scripts/pipeline/set_fc_reference.py`
   reads `calcs/dft/ferrocene` and **patches the shared `config/project.json`** (`reference.fc_abs_computed_V`).
   It has its own `fcntl` lock, but the real problem is upstream: each candidate's `native_thermal` runs
   `--only ferrocene,<gid>`, i.e. it **recomputes/writes the shared `calcs/dft/ferrocene/` namespace**, which the
   per-gid `flock` in motifs.py (L262, locks only candidate gids) does NOT cover. In a bag-of-tasks this is a
   cross-candidate write race on ferrocene. Mitigation: (i) make `fc_sync` a **single shared DAG node** depending
   on all `native_thermal:*` (as in §a), and (ii) hoist ferrocene to a **one-time prerequisite node**
   `native_thermal:ferrocene` that every candidate's `native_thermal:<id>` depends on OR that holds an exclusive
   `logs/pipeline_locks/ferrocene.lock`, so ferrocene is computed exactly once and never concurrently. This both
   removes the race and eliminates the redundant per-candidate ferrocene recompute.
3. **One task per GPU / gate-only GPU access**: every GPU node keeps `gpu_reserve run <gpus[slot]> -- taskset …`;
   `workers=len(gpus)` and the static slot→gpu map guarantee no two nodes target one GPU. Never set
   `CUDA_VISIBLE_DEVICES` directly (unchanged). Thread caps (`OMP/MKL/OPENBLAS/NUMEXPR/XTB`) stay pinned to
   `--threads` via `env` (L253-254) per worker block.
4. **Per-gid exclusive ownership** (L259-271) is retained; add the ferrocene lock from (c2).
5. **Certificate-locked kernels are untouched**: the native qRRHO/thermo kernel (`thermochemistry/native.py`,
   `NATIVE_THERMAL` with `validated_pyscf_version`/`validated_gpu4pyscf_version`, checked in `execution/gpu.py`
   `_probe_and_exec` L36-40), the projected-thermal certificate, and `config/native*.json` /
   `projected_thermal_validation.json` are read-only here. This refactor changes **scheduling and reporting
   only** — no change to argv content produced by `commands()`/`state_jobs()`, so every scientific command is
   byte-identical to today.
6. **Per-state isolation is what makes (a) safe**: `calcs/dft/<id>/<state>/` is already independent (documented at
   `state_jobs` L90-91), so the three dft states and distinct candidates never collide — the only shared writers
   are ferrocene/project.json, handled in (c2).
7. **Manifest invariants**: keep writing `run_manifest.json` incrementally and keep the terminal audit/verify
   gates (L319-363) so a failure still records `failed_stage`, `inspect_states`, and never sets
   `scientific_accuracy_validated=True` implicitly.

---

## Speedup estimate (GPU-utilization terms), 91 candidates

Rough per-candidate GPU-work: uma ~2.5 + dft(neu~18 + ox1~33 + ox2~18 = ~69) + sp ~0.25 + native ~5 + reorg ~6
≈ **~80 GPU-min/candidate** → ~91×80 ≈ **~7,280 GPU-min ≈ 121 GPU-h** of real work (native_thermal time is an
estimate — unmeasured; treat totals as order-of-magnitude).

- **Today, 3 GPUs:** only the dft fan uses 3 GPUs and it is ox1-gated (2 GPUs idle ~15 min/candidate); uma/sp/
  native/reorg/tables use 1 GPU while 2 idle. Wall ≈ ~42-45 min/candidate × 91 + 46-min tables ≈ **~64 h**,
  GPU utilization ≈ **~63 %** (121 GPU-h work / ~193 GPU-h delivered).
- **Bag-of-tasks, 3 GPUs:** 91 candidates × ~5 states-worth of GPU nodes ≫ 3 slots, so GPUs stay saturated;
  util → **~90 %**. Wall ≈ 121 GPU-h / (3 × 0.9) ≈ **~45 h**... corrected: ≈ **~42 h** (tables moved off path).
  **Speedup ≈ 1.5×** vs today, mostly from killing ox1-gating idle + single-GPU-stage idle + the tables tail.
- **6 GPUs:** today gains ~nothing (dft only has 3 states; GPUs 4-6 always idle) → still ~64 h.
  Bag-of-tasks on 6 GPUs ≈ 121 GPU-h / (6 × 0.9) ≈ **~22-24 h** → **~2.7-2.9× vs today's 3-GPU baseline**,
  ~1.9× over the improved 3-GPU case (near-linear: plenty of independent units).
  Caveat: `execute_graph` is **single-node**. 6 GPUs usually = 2 nodes × 3. Simplest: run **one `execute_graph`
  controller per node over a disjoint candidate half** (static split, like `launch.py` shards) — near-linear,
  each node ~45 candidates ≈ ~22 h, and the per-gid/ferrocene locks keep the two controllers from colliding on
  shared NFS (ferrocene must be computed by exactly one — assign it to one controller). A true cross-node dynamic
  GPU queue would mean a new claim-based layer (analogous to `chunked.py` but GPU-aware, going through the fleet
  dispatcher); defer unless the static split's imbalance proves costly.
- **Tables:** `--workers 8` turns the 46-min serial tail into ~8-12 min AND it is off the critical path, so it no
  longer blocks GPUs or risks discarding a successful QM campaign.

---

## Biggest risk
The **shared ferrocene / `config/project.json` state** (constraint c2). Bag-of-tasks makes many candidates'
`native_thermal`/`fc_reference_sync` want to touch the one `calcs/dft/ferrocene/` namespace and the one project
config concurrently — a write race the current per-gid lock does NOT cover and the serial loop accidentally hid.
It must be serialized (one-time ferrocene prerequisite node + a `ferrocene.lock`, and a single `fc_sync` node
after all native_thermal) before any parallel run, or provenance/Fc-reference corruption is possible.
