# Historical fleet sources

These `.sh.txt` files preserve the launchers reviewed on 2026-10-18. They are
reference text, not supported execution paths. Several use direct CUDA visibility,
incomplete CPU/RAM preflight, stale host policy or queue completion markers without
scientific verification.

Current paths:

- `python -m redox.workflows.resources` discovers availability through the system gate.
- `python -m redox.workflows.launch` plans DFT/UMA shards across suitable local GPUs.
- `python -m redox.workflows.gpu` preflights and reserves one full molecular QM command.
- `scripts/pipeline/run_motifs.py` orchestrates registered candidates with per-stage
  reservations and explicit scientific checks.

The old DFT/UMA wrappers now delegate to these implementations. Old claim-based
D3TaLES, reorganization, OROP and SCF/SP shell entry points stop with a migration
message. To run one of those validation calculations, use its Python validation
module through the generic GPU workflow after selecting inputs and reviewing the
current resumability and verification contract. Do not revive the old shell queues.
