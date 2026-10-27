# Reconciliation of twelve unstarted geometry intents

The [audit](audit.json) records process checks across all nine currently permitted
hosts, original dispatcher PID/start identities, the reservation-status hash and
an absence of scientific attempt directories and raw outputs. All twelve held
intents were verified unstarted. Each original intent is preserved beside the
audit before its replacement is marked `not_started`.

The live geometry controllers can now return these exact inputs to their normal
pending queue. No controller or scientific process was signaled, no completed
result was replaced, and reconciliation itself submitted no scientific jobs.
Unknown hosts, matching workers, live original owners, existing reservations or
newly appearing attempt files cause the reconciliation to stop without reopening
the affected case. Ten focused checks cover these guards and the inventory
adapter; see [the test report](tests.xml).

Live queue status remains in the original batches linked by the
[recovery audit](../geometry_queue_resumption/preflight_audit.json). This repairs
execution bookkeeping and does not qualify energies, free energies or ensembles.
