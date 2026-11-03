# Checked T4 chloride geometry resource class

The pilot now passes [independent geometry and resource review](independent_review.json):
102 actual SCF frames, a final maximum Cartesian gradient of 3.80282e-7 Ha/Bohr,
and 9,762 MiB peak GPU usage. The declared strict optimizer and fresh final
stationarity/spin/graph checks pass. This qualifies the exact geometry resource
class only; it supplies no Hessian minimum, native G or association energy.

The [current record](current.json) links all three released signed repairs.
Their [actual starts](followups/actual_starts.json) verify exclusive GPUs 0/1/2
on lambda5, six disjoint scientific CPU bindings, all six thread caps of two,
actual SCF entry and checked T4 FP64 runtimes. Per-node preflight, CPU/RAM and
aggregate durable-storage admission remain recorded. Each branch still needs
its own closed independent geometry review. These are the previously prepared
N0_below plus and N0_outer minus/plus seeds; the original failed minima and
prepared namespaces remain unchanged.

The declared T4 envelope is 12,288 MiB GPU, 32 GiB host RSS, two disjoint CPUs,
2 GiB durable output and four hours, with the existing 32-GiB free-RAM floor.
The [actual start](actual_start_verification.json) verifies one GPU resident,
worker PID/start, all six thread caps of two and separate controller CPU.
[Node resources](node_resources.json), [node import/version probe](node_probe.json)
and [aggregate storage admission](storage_admission.json) are retained.

The complete pilot review checked the exact final graph, strict fresh
stationarity, actual final energy and every SCF setting/frame/spin, measured
resource use and T4 architecture before releasing matching signed repairs.
The three new starts remain unqualified until their own checks pass. Hessians,
native G, association and general molecular workloads retain their separately
validated 32-GB-class gate.

Use `python -m redox.workflows.counterion_geometry_resources --help` for prepare
and launch and `python scripts/validation/review_counterion_geometry.py --help`
for the independent closed review. Frozen worker/profile sources and complete
per-evaluation checkpoints are on durable shared storage at the raw root in
[current.json](current.json). The [31 focused checks](scope_checks.xml) include
both T4 roles, unchanged full-workflow rejection, exact ionic connectivity and
runtime paths. This resource pilot supplies no material or reaction-yield claim.
