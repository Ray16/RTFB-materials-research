# Bounded remote geometry placement

Repeated remote-placement RPCs timed out without scientific attempts. Their
cause was not established. The new adapter uses a nonblocking host CPU-allocation
lock, records remote placement phases, and stops assigning further inputs to a
host after an unverified RPC response. The original per-host inventory, CPU/RAM
checks, GPU suitability, atomic reservation and scientific worker are retained.
A busy host lock returns a resource wait without consuming an input attempt.

The previous nine exact CPU controllers were temporarily stopped and replaced;
scientific workers were not signaled. Three unstarted intents passed the
all-host process/reservation checks in the [audit](audit.json), and their
original bytes are preserved beside it. The [verification](verification.json)
records nine live replacement controllers with disjoint CPU allocations and
one-thread native pools. Together with the completed tenth batch, the snapshot
has 77 checked geometries and 113 pending inputs, with no unresolved ownership
holds. No new scientific start or remote-placement receipt was established at
that snapshot.

Thirteen focused tests passed, including nonblocking lock contention, cleanup
on error, exact queue adaptation, fail-closed inventory and preservation of live
owners and existing scientific attempts. Commands, prior owner identities,
replacement launches and controller logs remain under
`calcs/validation/2026-10-26_geometry_bounded_recovery/`. Each controller continues
writing live status in its original batch; remote phases appear in each case's
`placement_checks/` directory when placement is attempted. This is execution
recovery, not an energy or ensemble result.

The [runtime placement audit](runtime_resource_waits.json) now records 72 actual returned resource waits. They returned in 0.010–10.674 seconds without scientific attempts: host allocation contention and a changed reservation snapshot remain normal queue waits. No new ownership hold is inferred from these returned calls.
