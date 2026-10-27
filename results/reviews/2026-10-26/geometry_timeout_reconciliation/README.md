# Follow-up reconciliation of three timed-out placements

Three remote placement requests timed out without a launch record, scientific
attempt directory or raw calculation output. Their exact CPU controllers were
verified by command, PID/start identity and process group, then stopped. Process
checks across all nine permitted hosts and the reservation gate established that
no corresponding worker or original-owner reservation remained. The original
intents are archived beside the [reconciliation audit](audit.json).

The same three controllers were restarted with explicit resume and distinct CPU
allocations. The [verification](verification.json) records all nine live geometry
controllers, their disjoint CPU sets and native thread limits. Across the ten
expanded batches, 77 geometries have passed independent review and 113 inputs
are pending; no ownership hold remained at that snapshot. The tenth dispatcher
is terminal because all nineteen of its geometries are already checked.

Controller logs and before/relaunch records remain under
`calcs/validation/2026-10-26_geometry_timeout_recovery/`. The separate independent
reviewer remains live. Geometry completion does not qualify electronic energies,
native free energies or ensemble convergence. The cause of the transient remote
placement timeout was not established; the timed-out requests and their errors
remain preserved in their original case directories.
