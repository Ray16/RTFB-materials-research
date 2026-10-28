# Improvement backlog (quality + parallelization)

From a three-part read-only review of the refactored tree (QM compute, orchestration/
parallelism, science/data). All findings were verified against the code. The suite stays
green; these are efficiency / parallel-tail / maintainability items, not correctness bugs.

## Implemented

- **Manifest read caching** — `core/common.py:read_manifest` now caches by (path, mtime,
  size), eliminating repeated NFS CSV re-parses in the hot per-state read path. Invalidates
  on any rewrite (test: `test_read_manifest_cache_never_serves_stale_rows`).
- **Bounded-memory evidence hashing** — `validation/spin_ensemble.py:digest` uses
  `core.hashing.sha256_file` (streamed) instead of `Path.read_bytes()`, so hashing a large
  `hessian_cartesian.npy` no longer loads it whole into RAM. Identical digests.
- **Drop redundant read** — `properties/potentials.py` reuses already-loaded `rO`/`rR`
  records instead of a third `read_result` call in the INCOMPLETE status string.
- **Unified atomic writes** — `validation/spins.py:atomic_json` delegates to
  `core.io.atomic_write_json` (fsync + `allow_nan=False`), matching the writer already used
  elsewhere in that module; no more silent invalid-`NaN` JSON.

## Deferred (higher value, higher risk — need GPU validation or design decisions)

- **[HIGH] Dynamic claim-queue sharding for the DFT/UMA launcher.**
  `workflows/launch.py` + `qm/dft.py` split the manifest by a static `k % n` stride. QM cost
  per molecule varies by orders of magnitude, so GPUs that finish their stride sit idle while
  the slowest shard runs. Every other stage (`chunked`, `benchmark_native`, `revalidate_sp`)
  already uses an atomic claim queue; the launcher should too, so a freed GPU immediately
  pulls the next unclaimed state. Biggest wasted-GPU-time item. Deferred: changes live launch
  semantics and is not unit-testable without the cluster.

- **[MED] Default CPU-thread caps in `dft.main`/`sp.main`.** Unlike the native/projected entry
  points, these don't cap pyscf/OMP/MKL threads by default. The proper launcher pins per-worker
  affinity via `taskset`, so this only bites when `dft`/`sp` run standalone; a wrong default
  could also hurt legitimate single-node runs. Decide the intended default before changing.

- **[MED] Parallelize per-candidate aggregate/RMSD.** `properties/reorg.py:aggregate()` and
  `potentials.main` run independent, pure-read RDKit automorphism-RMSD work serially. A
  `ProcessPool` with `OMP/MKL/OPENBLAS=1` per worker (shared-node courtesy) would speed this
  up with no change to results.

- **[MED] Replace monkeypatch-driven dispatch with explicit APIs.** `uma_current.run_dft/
  run_uma` swap `dft.*`/`uma.*` module globals + `sys.argv`; `hessian_resume` patches a
  gpu4pyscf class method; `native_thermal` overwrites `get_avail_mem` (never restored, budget
  keyed on caller-frame-name introspection); `polymer_campaign` patches `dft` internals;
  `launch_current` splices argv positionally. Safe across processes today but fragile. Give
  `dft`/`uma` a `run_batch(..., geom_resolver=, mult_resolver=, out_root=)` API.

- **[MED] Startup-ramp cleanups.** `gpu.py` holds the host-wide flock across nvidia-smi +
  `process_iter` + a 60 s `gpu_reserve list`; preflight (nvidia-smi + full import matrix)
  repeats per worker; `motifs.py` does ~5 reserve→release cycles per molecule instead of one
  long-lived lease. All ramp-up latency, none risk gate bypass.

- **[LOW] Consistency cleanups.** De-duplicate `to_float`/`_rmse` (import from
  `core.common`); consolidate the two equivalent two-spin Boltzmann implementations
  (`spins.two_state_ensemble` vs `spin_ensemble.equilibrium`); narrow bare
  `except Exception: pass` in `reorg.py`; wrap bare `open()` in reports in `with`;
  cache `core.common.load_config` (mtime-keyed) to stop re-exec'ing config `.py` per call.

- **[LOW] Layering.** `qm/hessian_resume` and `qm/projected_thermal` import from
  `workflows/*`; move shared digest/verify/allocation helpers into `core` so QM doesn't depend
  on orchestration.
