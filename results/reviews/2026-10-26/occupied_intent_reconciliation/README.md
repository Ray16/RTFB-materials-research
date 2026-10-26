# Occupied-chain submission-intent reconciliation

[Checked snapshot](summary.json) · [Ownership audit](audit.json) · [Tests](tests.xml) ·
[Prior controller startup](../occupied_followup_resumption/README.md)

Twelve prior electronic submissions were marked ownership-unverified after a
multi-host inventory timeout. The original controller PID/start identities were
absent, and checks on its four original hosts found no process referencing any
of the twelve exact cases. None had a launch record, attempt directory, worker
guard, scientific stage or raw calculation directory. Reservation status was
also checked. If any of these checks had been uncertain, reconciliation would
have stopped without reopening the affected entries.

The three replacement CPU controllers were paused during the checks. The original
intent files were copied and hashed here and on shared storage; the live intent
markers were then marked `not_started`, and the controllers were resumed. Their
next dispatch snapshots recognize all twelve as pending. At the recorded check,
all four controllers were live, no new scientific worker was running, and no
suitable GPU was reported free. Counts may change; the live batch dispatch and
inventory files, not this frozen snapshot, are authoritative for later progress.

This is ownership and scheduling recovery only. It establishes no new quantum
result, conformer convergence or candidate eligibility. Eleven focused tests
cover attempt-artifact rejection, remote-query handling and inventory behavior.
