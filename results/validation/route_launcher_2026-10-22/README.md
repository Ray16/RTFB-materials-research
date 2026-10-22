# Shared CPU launcher handoff

The [checked record](summary.json) verifies an actual launch through
`redox.workflows.cpu_chunks`. The launcher now dispatches the registered
`exact_routes` stage to its route worker; previously its command always invoked
the descriptor worker. Unknown and quantum stages fail before submission.

Two workers on lambda2 used disjoint one-core allocations with all six native
thread pools capped at one. They completed 31 and 22 chunks and checked all 837
exact designs. All remain held for missing exact-route evidence. No quantum jobs
were submitted. The [launch record](launch/launch.json), node preflight and
[source snapshot](source_snapshot/) preserve the actual execution epoch.

Raw inputs, checked outputs and receipts remain at
`calcs/scale_validation/2026-10-22_route_launcher/`. Earlier route-stage results
and their frozen sources remain unchanged; this is a new execution epoch.
