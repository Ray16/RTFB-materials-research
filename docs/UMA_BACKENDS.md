# Current UMA execution

The installed `redox` command now uses the current entry points. New UMA work
prefers Nitro on supported GPUs and selects the retained FairChem 2.21/Torch
2.8 engine explicitly on V100. The original replay modules and cached results
remain available. No dependency versions changed when the public entry point
was installed.

```bash
redox uma --show-policy
redox launch --module uma --out logs/uma_campaign
redox launch --module uma --out logs/uma_campaign --execute
redox launch --module dft --out logs/dft_campaign --execute
```

The launcher keeps the reservation gate, disjoint CPU sets and declared shards.
Perform the existing precursor/route review before launching a new full-property
campaign. `--require-nitro` refuses the V100 fallback. Actual Nitro requires
compute capability 8.x–12.x and its separately pinned interpreter; the dispatcher
re-executes that interpreter in the reserved process. The CUDA13 environment
cannot execute on V100. See the [policy](../config/uma_execution.json).

Composition, explicit charge and spin form the calculator state. A new fixed-state
predictor is created on a state change. ASE cache checks include charge/spin even
at identical coordinates, and an explicit reset also resets the inner calculator.
Per-state evaluation checkpoints contain coordinates, energy, forces and state;
they remain unqualified until optimization and graph/stereo checks pass.

Nitro results live under `calcs/uma_profiles/nitro/`. Model identity includes the
official current loader and execution policy. Historical results under `calcs/uma/`
retain their original identity. The changed composition-averaging semantics are
documented in the [V100 diagnosis](../results/validation/uma_v100_2026-10-24/README.md);
no diagnostic model-ID override is used. The DFT entry point prefers a complete,
graph-checked current warm start, honors an independently checked DFT spin override,
and retains the existing source when no new warm start exists. A partial or stale
Nitro record causes an explicit error.

These are gas-phase starting-geometry and diagnostic engines. Solution energies
and native thermochemistry use the existing DFT/SMD and DFT-Hessian protocols.
The [current validation](../results/validation/uma_current_2026-10-24/README.md)
passes the V100 dispatcher comparison on ten states and forty geometries. The
fixed-device official FP32/Nitro comparison is queued for a free H100; no redox
Nitro speedup or full optimizer migration is inferred from installation or the
V100 fallback test.
