# Efficiency / method-improvement survey (read-only)

Distinct from the already-planned scheduler bag-of-tasks, ferrocene fix,
incremental review_tables, and the λ_i self-consistent-surface fix.

## Executive summary — top 3
1. **Gate native DFT qRRHO to Pareto leads only.** Effort S. Native thermal
   (analytic DFT Hessian + ≥4 finite-difference SCFs + 3 xtb thermo passes per
   state) is the most expensive confirm step and today runs for every candidate.
   Driving it from a leads list (σ-aware Pareto front per pool, top-N per family)
   = ~6–9× fewer native Hessians. NO number change; certificate-safe
   (selection-layer only, native.py untouched).
2. **UMA physics triage over the large registry before DFT selection.** Effort M.
   discovery.py picks DFT candidates from cheminformatic descriptors only (MW, SA
   score, family rules) — zero energetics. A UMA gas-phase ΔE (approx IP/EA)
   prescreen (~27–290 s/state vs ~40 min DFT) can rank the 912-graph expansion
   and forward only the in-window band → 2–5× fewer DFT candidates on the live
   expand campaign. Changes *which* candidates enter DFT (selection), not the DFT
   numbers → needs a one-time UMA-vs-DFT calibration on the existing 91 to set
   the window + margin. Coverage risk → generous margin, prioritizer not filter.
3. **Multi-node fleet scale-out.** Effort M. launch.py already fans static
   `--shard N:i` across one host's GPUs; the bag-of-tasks plan is single-node.
   Run across free lambda0–6 GPUs via lambda-fleet → near-linear 2–4× over one
   node. NO number change. Risk: the shared ferrocene/project.json NFS write race
   is WORSE across nodes (one node computes ferrocene once under a shared lock).

## OFF-LIMITS (certificate-locked — verified by reading source)
- Named in SHARED_CONTEXT: qm/qrrho.py, qm/xtb_kernel.py, qm/xtb_runtime.py,
  quality/sampling_snapshot.py, hash-locked example payloads.
- **ALSO effectively certificate-locked:** thermochemistry/native.py carries
  PRODUCER_SOURCE_SHA256 (self-hash at import) and validates
  config/native_thermal_validation.json against a content sha256 + pinned
  pyscf 2.14.0 / gpu4pyscf 1.8.1. thermochemistry/projected.py and
  config/native*.json / projected_thermal_validation.json are the same
  certificate. Native SCF settings (grid_level=3, conv_tol=1e-12,
  df_batch_cap_GiB=6, finite-difference steps) are OFF-LIMITS.
- backends/screening_barrier.py is the "corrected" GPU grid-screening kernel
  wired to the native certificate; sensitive — do not touch for SCF speed.
- NOT locked but changing them changes PRODUCTION ENERGIES → requalify against
  the benchmark: backends/dft.py SCF guesses/grids/DF/opt protocol; sp.py adoption.

## Full opportunity list (value/risk ranked)
1. **E — native-leads gating.** S. No number change. Certificate-safe. Biggest
   confirm-cost cut (~6–9×). native.py's fcntl lock/circuit-breaker unchanged.
2. **A — UMA registry triage.** M. Selection-only (+calibration). Biggest
   absolute saver for the expand campaign. Never let a UMA ΔE appear as a
   reported potential (gas-phase, no certificate).
3. **F — multi-node fleet.** M (S for static split; L for dynamic cross-node GPU
   claim queue). No number change. Ferrocene-race caveat amplified; node-local
   env drift (launch.py/lambda-fleet preflight).
4. **B — UMA-supplied initial Hessian for DFT opt.** M. Neutral for rigid
   aromatic cores (shorter path, same converged geometry); floppy species can
   land in a different minimum → requalify-lite, opt-in + spot-check.
5. **C — closed-shell SCF warm-start.** M, LOW PRIORITY. Closed-shell SP seeded
   from opt-level density = identical energy, but SP is already cheap (~15 s;
   sp.py mostly ADOPTS). Open-shell 3-guess sweep (minao/atom/huckel) is a
   deliberate accuracy safeguard — do NOT touch (documented 2.5 eV / 0.11 eV
   wrong-guess failures).
6. **D — conformer/UMA-scan pruning.** M, LOWER PRIORITY. DFT already defaults
   nconf=1; real consumer is uma.py relax_state_scan (up to ~10 confs × mults ×
   states). Pruning risks spin-decision regressions → requalify-lite. Tune only
   if UMA wall becomes a bottleneck.

## Already-optimal (no action)
- RI-J density fitting ON for opt + SP (~3–5×, sub-meV cancellation).
- sp.py adopts protocol-matching energies without recompute (provenance guard
  correct — do not loosen).
- COSX/RI-JK for hybrid-exchange SP would change energies → not recommended (SP
  already cheap via adoption).
