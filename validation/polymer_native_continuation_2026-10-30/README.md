# Twelve occupied-state native continuations

[Summary](summary.json) · [Actual worker inspection](twelve_actual_workers.json) ·
[Exact inputs](inputs.json) · [Launch records](launches.json) · [Tests](tests.xml)

Twelve previously interrupted site-2 native states have actual scientific workers
on separate V100-32GB GPUs: three on lambda1, five on lambda2 and four on lambda4.
The inspection checks PID/start identity, direct-parent GPU lease ownership,
actual visible GPU, disjoint two-CPU allocations and six thread-pool limits of
two. Node capacity and remaining RAM promises are checked before placement; each
worker retains its 80 GiB host bound and 32 GiB free-memory floor. An available
GPU is insufficient when another worker already has a RAM promise.

The archived source, exact plan bytes, coordinates and scientific settings are
unchanged. A read-only filesystem namespace presents the archive at its original
checkout path. Only the fresh campaign directory, fresh raw-output parent,
temporary directory and authoritative shared reservation directory are writable.
The reservation tool runs inside that namespace so the original direct-parent
lease check remains intact. It is not bypassed or weakened. The archived worker
still validates its runtime, derivative, minimum, gradient and resource checks.
The new namespace source is captured with each launch; original attempts remain
at their archived locations. New attempts have distinct physical output paths
listed in [inputs.json](inputs.json).

Original-source/input hashes and CPU dependencies are checked inside each node's
namespace before submission. The runtime verifies the original PySCF 2.14.0,
GPU4PySCF 1.8.1 and cuTENSOR 2.2 identities and its FP64 contractions under a live
GPU reservation. The retained model is the separate smooth-surface acetonitrile
control; its results cannot silently replace the primary-model PMDI cycle.

Startup records preserve a rejected gate-outside-namespace experiment and an
occupied-device refusal. An initial inference from an exited wrapper PID was
corrected after checking the surviving gate and actual predictor. No restart is
justified by wrapper disappearance alone. H100 preflight initially found no
namespace utility; a checksum-identical user-space copy was staged and tested.
Subsequently foreign H100 residency caused seven placement refusals and one
owned-process interruption before native science. The affected controllers were
confirmed stopped before fresh V100 continuations. Foreign processes were not
signaled. These failures remain distinct from quantum nonconvergence.

The unchanged archived independent reviewer is invoked only after a native state
finishes and its exact owner and residency guard are closed. The active observer
and input list are bound in [active_observer_launch.json](active_observer_launch.json)
and [active_observer_inputs.json](active_observer_inputs.json). Its live report is
`results/validation/polymer_namespace_2026-10-30/active_independent_reviews/status.json`.
Twenty-four focused software checks pass. **No completed native free energy is
claimed by this startup milestone.** Failed or unfinished attempts remain outside
ensemble or material admission.

Preparation and submission use
`python -m redox.workflows.archived_native_namespace prepare`, then the node-local
`probe` and `launch` actions. The launch CLI must remain alive until its namespace
child returns; fleet submissions detach that persistent controller and preserve
its output. The caller must verify prior owner/predictor termination and the
current system host policy before submission. The versioned fleet and observer
entry points are `scripts/pipeline/launch_archived_native_followup.py` and
`scripts/pipeline/observe_archived_native_namespace.py`. They use only the declared
exact cases and do not enumerate new chemistry. More upstream and occupied-graph
work is still needed to close the polymer sampling review.
