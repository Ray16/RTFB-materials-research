# Prospective T4 chloride geometry resource pilot

The [current record](current.json) documents one actual exclusive optimization
start on lambda5 GPU 0, after the independently checked [fixed-gradient
witness](../counterion_T4_resource_2026-11-02/fixed_gradient_review.json).
This is the previously prepared N0_below minus signed-mode repair: exact
42-atom/412-AO oxidized viologen core plus chloride, total +1 singlet.
Its primary r2SCAN-D4/def2-SVP legacy-SMD MeCN method, corrected integration,
strict optimizer convergence and graph/stationarity checks are unchanged.
The original prepared namespace and failed minimum remain preserved.

The declared T4 envelope is 12,288 MiB GPU, 32 GiB host RSS, two disjoint CPUs,
2 GiB durable output and four hours, with the existing 32-GiB free-RAM floor.
The [actual start](actual_start_verification.json) verifies one GPU resident,
worker PID/start, all six thread caps of two and separate controller CPU.
[Node resources](node_resources.json), [node import/version probe](node_probe.json)
and [aggregate storage admission](storage_admission.json) are retained.

Optimization is running and remains unqualified. A closed, independent review
must check the exact final graph, strict fresh stationarity, actual final energy
and all SCF settings/frames/spins, measured resource use and T4 architecture.
Only a passing full geometry review may release the other prepared signed
repairs. No extra repair is released now. Hessians, native G, association and
general molecular workloads retain their separately validated 32-GB-class gate.

Use `python -m redox.workflows.counterion_geometry_resources --help` for prepare
and launch and `python scripts/validation/review_counterion_geometry.py --help`
for the independent closed review. Frozen worker/profile sources and complete
per-evaluation checkpoints are on durable shared storage at the raw root in
[current.json](current.json). The [31 focused checks](scope_checks.xml) include
both T4 roles, unchanged full-workflow rejection, exact ionic connectivity and
runtime paths. This resource pilot supplies no material or reaction-yield claim.
