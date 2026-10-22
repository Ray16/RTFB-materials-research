# Source-preserving last-spin continuation

The first actual alternate-triplet run exhausted a V100 32-GB device after its
geometry optimizer converged. Its original failure, trajectory and real optimizer
return remain in `../spin_completion_2026-10-21/`. Neither a callback nor a saved
XYZ by itself qualifies an optimized geometry.

The [plan](plan.json) pins the actual converged optimizer return, its matching
final SCF evaluation, the original chosen physical state and numerical producers.
The new driver reuses only this exact helper geometry and performs the ordinary
active energy call in a fresh process, without retaining optimization SCF objects.
It records that provenance explicitly. The normal source-addressed spin comparison,
topology, SCF, active SP and any required native-minimum checks remain in force.
This is not quantum restart-equivalence evidence or a new synthesis candidate.

The finite CPU coordinator searches all policy-allowed nodes and requests one
free device with at least 80 GiB for this one remaining task. An occupied device
is unavailable. It uses the standard GPU wrapper for dependency, CPU, memory and
reservation checks. Contention does not consume an execution attempt; a real
worker failure ends this continuation for review. The resource window is 24 hours.

Read [coordinator.json](coordinator.json) and, once present, `status.json` for actual
progress. Prepared input and a waiting coordinator do not count as a running GPU
calculation. A completed energy is not automatically spin admission. Three
negative/continuation tests pass in [tests.xml](tests.xml), checking convergence,
geometry identity, exact spin scope, real scoring and failure preservation.
