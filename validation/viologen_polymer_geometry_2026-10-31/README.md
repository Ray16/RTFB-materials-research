# Exact viologen occupancy geometry campaign

Ten resting-state ionic polymer models now have verified live starts in MeCN:
three single grafts and adjacent/one-spacer pairs in each of the all-same and
alternating backbone configurations. The [captured observation](startup_observation.json)
checks actual PID/start identity, frozen worker arguments, exclusive GPU residency,
disjoint two-CPU affinities, six native thread caps, the NumPy allocation profile
and actual SCF iteration progress. Eight V100 workers are on lambda1, one on
lambda2 and one on lambda4. Nine are newly submitted; the largest original
optimizer is retained in its original namespace.

This is an execution milestone. None of the ten structures has been accepted as
a minimum or converged conformer ensemble in this snapshot. No charged-polymer
free energy, occupancy nonadditivity, electron-transfer coupling or all-axis
material decision is reported. It advances M1 without closing M1 or M4.

Each graph has its prospective exact-route check and initial XYZ topology/stereo
audit. Coordinates represent the bare ionic polymer with omitted separated
spectator ions. PF6 dry-salt mass bookkeeping is separate. An unoccupied reactive
site retains benzyl chloride; terminal plain-styrene padding remains distinct.
The solution calculation uses the explicitly named SWIG/CDS alternate SMD model,
r2SCAN-D4/def2-SVP, the isolated grid-screening correction and the compact 28-GiB
execution profile already compared against a direct large-GPU initial-state
reference. It does not qualify the complete primary-SMD reaction cycle.

Every evaluated geometry, energy, gradient and SCF convergence flag is retained.
Rolling current and final SCF checkpoints bound each optimizer's durable output
to 256 MiB, with an eight-hour wall budget. Aggregate admission also reserves the
remaining output budgets of the three running full-Hessian controls and a 32-GiB
shared-storage margin. The legacy full native-thermochemistry storage gate is
unchanged. All suitable free GPUs needed by this roster are considered; no small
worker cap is imposed.

The new independent final-geometry reviewer binds optimizer-returned coordinates,
XYZ and the actual final SCF checkpoint; checks fresh per-atom force norms,
charge/spin, retained indexed covalent connectivity, stereochemistry and bond
lengths; and checks normal completion, resources and ownership. It deliberately
leaves Hessian minimum and ensemble/free-energy gates incomplete. The 62 focused
checks include corrupted state/coordinate/energy, norm, stereo, graph, ownership
and storage-admission cases. All ten selected initial geometries also pass the
new geometry auditor, which is distinct from reviewing an optimized structure.

The [prepared roster](prepared_panel.json), [dispatch record](dispatch.json),
[summary](summary.json), [test report](tests.xml) and producer snapshots are frozen
by [manifest](manifest.json). Live outputs remain under the raw root recorded in
the summary. These running-start counts describe the capture time only.

To refresh process observations, run `scripts/validation/observe_charged_geometry_panel.py`
with the raw manifest and a new output path. Review a normally closed optimizer
with `scripts/validation/review_charged_polymer_optimization.py --spec NAMESPACE
--out NEW_REVIEW.json`. Preserve failed/uncertain attempts; do not reclaim them by
age or infer scientific success from a launch receipt.
