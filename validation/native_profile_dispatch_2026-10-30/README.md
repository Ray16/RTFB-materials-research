# Matched resource profiles in the occupied-state workflow

The native dispatcher now uses an independently checked prospective memory
profile for matching pending tasks. Seven profiles were prepared from the
completed `all_same__occ_4__conf003` pilot; two additional calculations actually
started on V100-SXM2-32GB GPUs, one each on lambda4 and lambda1. Both have written
SCF checkpoints and base gradients. These starts are not completed native free
energies or a qualified occupied ensemble.

The new scheduler verifies the reference's native result and reconstructs its
resource qualification. Matching requires the same constitution, atom/basis
counts, charge, spin, native settings, solvent model and runtime. Each selected
task retains the original immutable 80-GiB plan and inner guard and adds the
prospectively declared 40-GiB outer guard. Existing jobs retain their original
reservations. Other task classes use the ordinary resource bound. This measured
neutral-singlet class does not qualify charged candidates or every oligomer.

Scheduling scans pending tasks for a fitting allocation, so an unrelated larger
task at the front cannot block a matching smaller task. The electronic-energy
dispatcher defers only when a real pending native task fits the available
resources. Shared GPU reservations, per-node dependency probes, two disjoint CPU
slots, six thread caps, host load/RAM checks, unused RAM promises and the
three-failure circuit breaker remain in place. No fixed GPU-worker cap is added.

Completed wrapped jobs follow the original numerical reviewer with its exact
execution-wrapper identity check, then the separate effective-resource review.
The completed pilot reconstructs through this dispatcher path, reproducing all
three rotor-cutoff free energies and 2,449 source bindings. Forty-one focused
tests cover immutable scheduling copies, source/reference mismatches, budget
headroom, head-of-queue blocking, exact priority ownership, changed inputs and
the shared-volume admission floor.

The shared output volume had approximately 91 GiB free during preflight. New
native and electronic starts now require at least 64 GiB free; running workers
are not interrupted by this admission check. This is a conservative shared-volume
headroom policy, not a measured per-job disk guarantee. The alternate mounted
volume has no writable rzhu directory at the inspected locations, and a durable
alternative has been requested. Existing source evidence and live outputs have
not been moved or removed.

Only the two owned CPU dispatchers were stopped for the update, with exact
PID/start checks and no direct subprocess at termination. Detached GPU workers
continued. The first stop pass terminated the energy controller and deferred
the native stop while a subprocess was active; the second pass confirmed both
CPU controllers stopped before editing. Old receipts and statuses are preserved.

Current runtime outputs are:

| Stage | Output |
| --- | --- |
| Native dispatch and checked results | `results/validation/native_conveyor_2026-10-30/dispatch/status.json` |
| Native preparation stream | `results/validation/native_conveyor_2026-10-30/preparation/status.json` |
| Electronic dispatch | `results/validation/energy_classes_2026-10-30/queue/status.json` |
| This update's owner/start evidence | `results/validation/native_profile_dispatch_2026-10-30/` |
| Raw native inputs/checkpoints | `/nfs/lambda_stor_01/homes/rzhu/redox_polymer_continuations/2026-10-30_native_conveyor/` |

[summary.json](summary.json) and the two dispatcher snapshots describe the
checkpoint, not later live execution. [actual_scientific_starts.json](actual_scientific_starts.json)
verifies exact PIDs, GPU residency, affinity and thread caps.
[scientific_checkpoint_progress.json](scientific_checkpoint_progress.json)
records the actual durable gradient/SCF output. The frozen paper's 75-state
snapshot remains unchanged; new occupied results enter the scientific analysis
after their physical, duplicate and ensemble checks. M4 remains open.
