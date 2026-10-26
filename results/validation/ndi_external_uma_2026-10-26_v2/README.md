# Isolated NDI external-control UMA campaign

[Frozen campaign plan](plan.json) · [Adapter tests and startup snapshot](../../reviews/2026-10-26/ndi_uma_adapter/README.md) ·
[Source and experimental limits](../ndi_external_challenge_2026-10-26/README.md)

The source-checked NDI control has an isolated three-state manifest and ten
neutral-graph starting geometries. `scripts/pipeline/uma_external_control.py`
stages those checked inputs and uses the existing state-aware UMA runner with a
process-local library override and separate output directories. The shared
candidate manifest and existing candidate caches are not changed. State charge
and multiplicity values are initial hypotheses, not completed state assignments.

The worker requires an active `gpu_reserve` lease, selects Nitro on supported
architectures and the previously validated legacy fallback on V100, records
claims and failures, and keeps the existing optimizer/checkpoint and geometry
checks. The scheduler must select the matching interpreter and library path
before reservation; the worker refuses an interpreter mismatch rather than
replacing a live residency guard. A direct invocation without a reservation was rejected. Mocked tests
checked manifest isolation, output routing, failure preservation and restoration
of the shared module state. They did not execute a model or validate the
external molecule's energies.

**No UMA or DFT state has started for this control.** The fixed-device
redox-specific FP32/Nitro comparison is running as a separate resource waiter;
its numerical qualification remains pending. The V100 fallback comparison
and the H100 Nitro pair have different scopes. Solution DFT, alternative-spin
checks, native thermochemistry and comparison with the peak-based experimental
proxies are still required. A [neutral-state resource-pilot coordinator](../../reviews/2026-10-26/nitro_pair_independent_observer/README.md)
is live and gated on independent Nitro parity. Fan-out to the other states and
the solution DFT handoff remain to be implemented after the pilot is reviewed.
