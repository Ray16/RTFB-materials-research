# Native conformer controller: per-host inventory recovery

[Checked snapshot](summary.json) · [Inventory snapshot](inventory_snapshot.json) ·
[Dispatcher snapshot](dispatch_snapshot.json) · [Focused tests](tests.xml) ·
[Live dispatcher](../../../polymer_models/campaigns/2026-10-26_native_final_state_witness_batch/dispatch.json) ·
[Live inventory](../../../polymer_models/campaigns/2026-10-26_native_final_state_witness_batch/resilient_inventory/status.json)

The previous controller exited when a multi-host availability query timed out.
Its completed results and prior dispatch record remain preserved. The recovered
controller uses a separate, source-bound adapter that queries allowed hosts
individually. A timeout or malformed response is recorded as an **unverified
host** and cannot supply a free GPU; successful host inventories remain usable.
The original atomic GPU reservation, node preflight, scientific worker and
physical gates are unchanged.

At the recorded snapshot the controller was live, pinned to one disjoint CPU,
with four native free energies checked and no scientific worker running. All
nine configured allowed hosts gave verified inventory responses and reported
no available GPUs at that time. Availability can change; use the live files and
actual ownership records rather than interpreting the snapshot as a current
resource allocation. The controller monitors for suitable free devices within
its finite 24-hour window. A controller start is not a scientific calculation
start or a converged ensemble.

The two focused tests exercise partial-host timeouts, correctly attributed free
slots and rejection of malformed, duplicate or misattributed results. No
scientific output was relabeled or moved. The separate experimental solvent
model remains distinct from the primary manuscript calculations.

Source: `scripts/pipeline/coordinate_native_inventory_resilient.py`.
Run log and launcher ownership:
`calcs/validation/2026-10-26_local_motif/native_resilient_20261026/`.
