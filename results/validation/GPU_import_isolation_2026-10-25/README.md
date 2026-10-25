# CPU dependency checks and GPU import residency

The [reserved import probe](reserved_probe.json) found zero resident contexts
before importing GPU4PySCF and one 100 MiB context afterward. No explicit
quantum calculation was run. Therefore importing this library is not a GPU-free
CPU dependency check. The project now uses
`scripts/fleet/preflight_cpu_dependencies.py`: CPU modules are imported, GPU
distributions are checked through metadata, and a process-wide import guard
blocks direct/transitive GPU module imports. Actual GPU imports and runtime
checks remain inside each reservation, before scientific execution.

[The node check](CPU_preflight_lambda2_final.json) reports deferred GPU imports
and an empty loaded-GPU-module list. [CPU controller transitions](CPU_controller_transition.json)
and [restarts](CPU_safe_restarts.json) preserve the earlier controllers and
change only their dependency-probe adapter. Scientific producer bytes and live
workers were not changed or signaled. The canonical fleet skill is unchanged.

Two occupied-geometry attempts recorded foreign resident PIDs and their own
guards stopped them. Those historical PIDs are gone, so attributing those two
events specifically to a dependency probe is not proven. Both source attempts
remain excluded and preserved. Separate same-model continuations use their
durable unqualified geometry frames and repeat all normal checks under
`results/polymer_models/campaigns/2026-10-25_smooth_occupied_resident_recovery/`.
