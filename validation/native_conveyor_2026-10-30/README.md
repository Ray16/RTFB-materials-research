# Checked energies to independently reviewed thermochemistry

This continuation addresses the missing occupied-conformer free energies under
reviewer item M4. It advances the declared phthalimide control cohort through
electronic review, native-input preparation, guarded GPU execution and independent
thermochemistry reconstruction. It does not complete the occupancy ensembles or
establish a material shortlist.

[summary.json](summary.json) freezes this milestone's observed counts and worker
checks; [manifest.json](manifest.json) binds the implementation and small evidence.
The changing runtime pointers below are authoritative for subsequent progress.

| Stage | Current output | Acceptance |
| --- | --- | --- |
| Electronic evaluation | `results/validation/energy_classes_2026-10-30/independent_reviews/status.json` | Checked exact graph, geometry, MeCN model, raw energies and source hashes |
| Native inputs | `results/validation/native_conveyor_2026-10-30/preparation/status.json` | Identical electronic parent and geometry; existing namespaces reused without duplicate jobs |
| Native dispatch and physical review | `results/validation/native_conveyor_2026-10-30/dispatch/status.json` | Original SCF, gradient, Hessian, frequency and partition checks, followed by source-bound independent reconstruction |
| Smaller-memory pilot | `results/validation/effective_native_memory_2026-10-30/independent_reviews_v2/status.json` | Original numerical review plus the stricter guard, Linux high-water RSS and exact execution-wrapper chain |

The preparation producer and dispatchers have bounded observation windows and
durable ownership receipts. They use the system GPU reservation policy, suitable
32-GB-class or larger devices, disjoint two-CPU bindings, six native thread caps,
host load/RAM checks, and unused RAM promises. The ordinary native dispatcher
retains the original 80-GiB guard and 32-GiB free-memory floor. A separate,
prospective 40-GiB pilot adds a stricter outer guard. That one tested run now
[passes both physical and resource reviews](../polymer_effective_checked_2026-10-30/README.md).
Its completed collector is `independent_reviews_v3/status.json` under the same
runtime pilot directory. Wider matched-class use still needs explicit prospective
profiles. No in-flight reservation is retroactively shrunk.

Ready native calculations receive admission priority over new energy-only work
when they fit a free device and its host. Running calculations remain untouched;
there is no fixed fleet-size cap. An uncertain submission retains its claim. An
explicit GPU-gate refusal can reopen only after both owners are proven gone and
no scientific outputs exist. Scientific failures stay preserved; three consecutive
scientific failures stop further submissions on that host. A stopped dispatcher
requires explicit resume and exact prior-owner verification.

One initial remote handoff failed after successful preflight. A five-host process
audit and absence of controller/scientific outputs established that no job had
started. The original claim and placement are retained under the runtime
`handoff_recovery/` history. The revised launcher preserves remote spawn stdout
and stderr and ensures its output directory exists on the target node. The
original stderr was unavailable, so a specific root cause is not asserted.

Raw calculations and checkpoints remain under
`/nfs/lambda_stor_01/homes/rzhu/redox_polymer_continuations/`, including
`2026-10-30_native_conveyor/` and `2026-10-30_effective_native/`. Original archived
inputs are mounted read-only at their recorded paths. Only derivative provenance
paths are relocated; frozen scientific inputs, equations and settings are retained.

[guard_cpu_controls.json](guard_cpu_controls.json) records an actual normal exit
and a caught transient allocation whose sampled RSS alone would have missed the
peak. [review_adapter_check.json](review_adapter_check.json) verifies the real
archived review function changes only its execution-program identity predicate
and checks visibility of the exact review source inside the pilot's namespace.
The runtime records are never relabeled to make them look like another program.

The updated numerical results will enter the article and SI after complete,
balanced occupied ensembles pass their physical and sampling checks. The current
paper's forty empty-input comparisons remain the latest integrated ensemble
results. Energy-only completions and successful launch receipts supply no missing
free energy.
