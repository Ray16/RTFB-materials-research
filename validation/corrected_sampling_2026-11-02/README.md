# Complete corrected polymer-sampling cohort

The [largest native pilot](double_native_review.json) now passes independent
full-Hessian, sampled-gradient, actual-wavefunction and nuclear-partition review.
Its 17 SCF frames are checked, its lowest internal mode is 6.4952 cm⁻¹ and its
GPU/host peaks are 24,062 MiB / 32.08 GiB. This releases matching native inputs
under the unchanged 31-GiB-GPU / 80-GiB-host bound. The
[new native collection](native_collection_04.json) contains eight checked inputs;
the remaining native and ensemble quantities are pending. Later preparation and
live scheduling counts are separate from this frozen collection.

The corrected native stage now has a maintained resumable coordinator:
`python scripts/fleet/run_corrected_native_queue.py --help`. It retains the full
240-input objective, collects closed native calculations independently and
holds wider preparation/dispatch until the largest native pilot passes. Two
native receipts have also passed the [new context/component reader](native_receipt_crosscheck.json).
All three rotor sums, the exact energy/geometry parent, thermal definition,
charge/spin, model, standard state and physical/resource QC remain attached to
each admitted value. No occupied ensemble is yet qualified.

The controller uses the existing reservation/preflight dispatcher, aggregate
RAM/storage budgets and append-only peer registry. It never restarts a failed
scientific worker. Its source-change guard was exercised during development;
the stopped controller and history are retained, and the updated coordinator
was explicitly resumed only after exact process closure. Raw controller records
are under
`/nfs/lambda_stor_01/homes/rzhu/redox_polymer_continuations/2026-11-02_corrected_sampling_native_followups/`.
The [focused checks](native_queue_checks_05.xml) pass 72 tests, including actual
CPU SCF/checkpoint controls and rejection of altered native values, context,
partial rosters, premature release, later residency interruptions and terminal
failures that must be reported rather than treated as active work.

The [source inventory](input_inventory_complete.json) accounts for all 240
declared phthalimide-control inputs: six occupancy graphs, two explicit backbone
configurations and twenty starting IDs per cell. All 240 actual XYZ inputs pass
full graph/hydrogen connectivity and stereochemistry checks during preparation.
Starting geometries comprise 122 physically replayed native geometries, 107
independently checked geometries, eight previously prepared scored geometries and
three separately checked electronic-reference geometries. No historical E or G
is imported. The first two incomplete inventories remain preserved.

The [prepared manifest](prepared_manifest.json) freezes a fresh corrected
gas/SMD electronic cohort at
`/nfs/lambda_stor_01/homes/rzhu/redox_polymer_continuations/2026-11-02_corrected_sampling_energy_v4/`.
The [sampling contract](sampling_contract.json) retains the original twenty-start
IDs, conditional configurations, 10-to-20 comparisons, three rotor settings,
duplicate sensitivities and physical derivative/minimum limits. These are
method controls, not synthesis candidates or measured occupancy probabilities.

The [current frozen complete-roster collection](electronic_collection_11.json) has
89 independently checked electronic inputs and no review errors. Every one of
the 240 slots remains present; pending E/G values are null. The current
[PSS-profile manifest](PSS_prepared_manifest.json) supersedes the earlier
resource preparations while preserving them. No new native G or occupied
ensemble cycle is supplied by this electronic collection.

The separate [native collection](native_collection_04.json) retains all 240
slots with eight checked G inputs. Its [controller snapshot](native_queue_status_03.json)
and [process observation](native_followup_observation_01.json) distinguish
active coordination, the live largest pilot and physical qualification.

## Method and storage

The existing ωB97M-V/def2-TZVPD gas and MeCN single-point equations are unchanged.
The MeCN model is the separately declared PySCF SWIG/CDS smooth-surface model,
with the run-local grid-screening correction. Native followups must supply their
own r2SCAN-D4/def2-SVP nuclear terms before this cohort can supply G. These hybrid
electronic energies must not be replaced by geometry-level r2SCAN energies or
mixed with historical/primary-model components.

An explicit occupied-orbital SCF format retains the complete Slater determinant,
all occupations, all orbital energies and solver settings. It omits virtual
coefficients from disk, supports a density restart and explicitly does not
support a full response/Hessian restart. Each iteration is checkpointed, each
kernel has its own directory and the original multiple-guess/fallback solver
retains a final determinant. No solver or physical setting is changed.

[Actual saved-determinant comparison](occupied_equivalence.json) covers a
1193-AO/235-occupied-orbital polymer control and a 394-AO exact molecular state.
Occupied coefficients, occupations, orbital energies and total energies remain
bitwise identical, and both density differences are zero. This CPU storage
comparison does not requalify source convergence or calculate new E/G.
Actual CPU SCF equality, density restart, interruption retention, orbital-metric
corruption rejection and multiple-guess restoration also pass tests.

## Starts, resource checks and pending results

### Qualified PSS profile and released electronic cohort

All three new energy pilots pass independent closed-execution, actual orbital/
geometry and gas/SMD energy review:
[empty](PSS_empty_pilot_review.json), [single](PSS_single_pilot_review.json) and
[double](PSS_double_pilot_review.json). Their effective host peaks are
72.070, 111.047 and 166.150 GiB, within the declared 160/160/240-GiB classes with
25% headroom. Peak GPU usage is 16,792/18,804/19,962 MiB. The
[cross-profile repeats](memory_profile_electronic_comparison.json) reproduce
the earlier empty/single gas/SMD energies within 4.37e-11 eV.

The [controlled fork experiment](fork_memory_control.json) distinguishes shared
pages from a real 96-MiB private allocation. Linux
[PSS accounts for each process's proportional share of shared pages](https://docs.kernel.org/filesystems/proc.html).
The new `corrected_hybrid_energy_pss` profile limits the maximum of the main
process RSS high-water mark and summed descendant PSS. It retains raw tree RSS,
child identities, the original main-process 160/240-GiB bounds and 32-GiB free-RAM
floor. A real single-graft run recorded 221 GiB summed RSS and approximately
111 GiB PSS. The older native/RSS profile is unchanged. Independent reviewers
accept PSS only when explicitly declared and reconstruct its effective peak;
missing/invalid metrics cannot pass.

The [recovery receipt](density_restart_source.json) proves the stopped double
pilot's exact ownership closure and molecular/orbital frame. Its unqualified
checkpoint supplied only a starting density. The new solution solve reconverged
in six iterations, followed by its own gas solve. The per-phase policy preserves
both the legacy nominal guess and the actual `checkpoint_density` guess, with
the exact input hash; no convergence or property qualification transfers.

After all three pilots passed, the bounded dispatcher released the full cohort.
Its first pass submitted ten more energy jobs, which now all pass independent
review. Subsequent passes retain atomic claims, exact stopped-owner checks and
all pending inputs. A ready smaller job can fill remaining host RAM when the
queue's first larger job does not fit; RAM promised to starting workers is
accounted separately. This improves utilization without changing resource bounds.

### Native thermochemistry

The [three native pilot namespaces](native_pilot_manifest.json) run the original
smooth-model native equations under corrected screening. Their electronic inputs
are independently checked; their own r2SCAN-D4/SVP SCFs, complete Hessians,
sampled derivatives, projected modes and nuclear partitions remain mandatory.
All three center SCFs and matched SCF/CPHF grids have been saved. The
[empty native input](empty_native_review.json) now passes full independent
reconstruction: 17 actual SCF frames, a 9.9656-cm⁻¹ minimum, sampled derivatives
and its own three nuclear partitions. The [site-2 single graft](single_native_review.json)
also passes its own 17-frame, derivative, mode and nuclear reconstruction, with
a minimum frequency of 11.0972 cm⁻¹. The largest native pilot, complete occupied
sampling and ensemble cycles remain pending.

The dimension witness separates a 119-atom corrected full Hessian that completed
before its producer's wall stop from its subsequently completed derivative/mode
validation. It does not invent a normal construction exit or a retained full-
construction GPU peak. New complete native pilots must pass their own 80-GiB/
31-GiB resource and physical gates before the largest pilot can release further
native inputs. The numerical reviewer reuses the original scientific expressions,
adapts only the historical guard fields and adds exact execution, actual orbital,
grid and original-equation response checks. Full virtual-orbital checkpoints are
retained for native recovery; raw Krylov vectors are not part of this bulk profile.

The [current combined suite](../ionic_cycle_completion_2026-11-02/current_combined_checks.xml) passes 165 tests. It
includes actual CPU DFT/Hessian equality, mismatch rejection for response grids,
real fork/private-memory accounting and density-restart reconvergence. A shared
receipt reader checks context and exact raw E values before collection or native
preparation, rejecting corrupted scalar caches and mismatched source states.

The latest [process observation](full_campaign_observation_02.json) is separate
from these accepted values. Quantum workers and budget observers retain their
own PID/start identities; observation timeouts cannot trigger recovery. The
collector retains its as-used source bytes and checks each cached E against its
exact bound raw scalar before reuse.

### Autonomous queue and exclusivity recovery

The bounded [queue observation](queue_observation_02.json) and mutable raw ledger
at `electronic_queue/electronic_current.json` retain all 240 slots and refill
eligible GPUs after review. The controller has an exact PID/start identity,
one CPU, a 24-hour window, source-change stop, prelaunch cooldown/attempt limits
and a three-failure host circuit breaker. It never restarts a failed science job
by age. Its storage-peer ledger is append-only and includes the ionic followups.

One producer returned normally but recorded a foreign GPU resident. Independent
closure review rejected it; its numerical result remains unqualified. The
[successor manifest](recovery_manifest_01.json) preserves that namespace and
replaces only the same exact input slot with a fresh density-initialized solve.
No source qualification transfers. The first recovery preparation's copied
bytecode directory was caught before any GPU execution; copying/hashing now
excludes bytecode directories and preserves that failed preparation separately.

A later initial-residency query exited before any scientific worker or SCF
started. [The exact guard, wrapper and controller were proven closed](guard_preflight_recovery_source.json);
the [next successor manifest](recovery_manifest_02.json) replaces only that
unchanged input with a fresh namespace. Existing guard records remain intact.
A query timeout cannot authorize a restart while any recorded process is live
or unverified; dedicated tests enforce this distinction.

[Independent CPU watchdogs](watchdog_observation_01.json) now supplement the
historical in-process monitor. An exact UID/PID-start/plan-argv check precedes
every owned-process signal; an ignored SIGTERM can escalate to a kill of that
same owned worker. Foreign processes are never signaled and claims/reservations
are never cleared. A recorded external interruption prevents normal scientific
admission even if a runtime suppresses the first termination signal. Actual
TERM-ignoring-process and wrong-identity tests pass.

### Preserved first-stage observations

The following v4 observations document the earlier held-pilot stage. The
qualified PSS profile and released cohort above supersede those holds.

[The frozen start observation](start_observation_01.json) confirms three live
energy workers and the existing larger viologen native recovery by exact
host/PID/start identity. Energy pilots occupy λ1 GPUs 0/1 and λ2 GPU 0; the lead
recovery retains λ13 GPU 1. Each worker has disjoint two-CPU affinity and all six
native caps equal two. Actual GPU residency and entry into all three pilot SCF
kernels are separately recorded. These observations are dated evidence, not
permanent assertions of live state or completed physical results.

The [empty pilot review](empty_pilot_review.json) now independently accepts both
gas and MeCN electronic energies, with one actual occupied determinant for each.
Its maximum coordinate error is 1.78e-15 Bohr and occupied-orthonormality error
4.67e-15. Peak GPU memory is 16,768 MiB and host RSS is 72.1 GiB. This qualifies
one electronic input, not its new native G or an occupied ensemble.

The [single-graft pilot review](single_pilot_review.json) also accepts both gas
and MeCN electronic energies, with independent actual-frame and occupied-orbital
checks. Its peak GPU memory is 18,856 MiB. Thus two of the three electronic
pilots are checked; no corresponding new native G is yet qualified.

The double-graft pilot stopped at its owned host-memory guard. Its main-process
RSS peak is 166.09 GiB but its summed process-tree peak is 332.03 GiB, above the
declared 240-GiB bound. Its last iteration checkpoint is retained unqualified.
Transient shared/forked memory is a hypothesis to test before recovery; the
record does not establish an out-of-memory event or permit ignoring the guard.
The larger lead recovery continues separately. The
[then-current observation](progress_observation_03.json) records exact ownership and
nine completed response blocks for that lead at its observation time. A complete
response is still followed by full Hessian, displaced-gradient and minimum QC.

Only the empty, single- and double-graft energy pilots may start until all three
have independent closed-execution, molecular-frame, occupied-wavefunction,
gas/SMD convergence and resource reviews. Followup release also requires 25%
host-RAM headroom within the 160/240-GiB profiles. The default native guard
retains its 80-GiB ceiling. A launch record supplies no physical qualification.
The remaining 237 energy inputs are held; all corresponding
new native/ensemble qualification remains pending.

Preflight caught an inventory hidden by the archived source overlay; the input
is now frozen outside that overlay. A separate guard-argument refusal exposed
the native 80-GiB profile mismatch before scientific execution. Both earlier
namespaces are preserved. An empty-pilot RAM refusal was retried only after the
exact controller was proven stopped and no scientific worker/guard had started.
The original placement and receipts remain intact. Shared dispatch now accounts
for RAM promised to not-yet-started workers, permits only proven prelaunch retries
and preserves foreign/uncertain claims. As-used dispatcher snapshots and their
explicit source relocations are retained alongside these records.

The [focused checks](final_checks.xml) pass 120 tests, including real CPU SCF and
same-process memory-guard execution, asynchronous RAM assignment and refusal of
foreign/unverified recovery. Full quantum and occupied-sampling completion,
counterion/speciation sensitivity and publication release remain outstanding.

## Maintained entry points

Use `python -m redox.workflows.sampling_energy --help` for preparation/dispatch;
`--pilots-only` restricts the initial stage, and full dispatch requires three
exact `--pilot-review` receipts. Use
`scripts/validation/review_corrected_sampling_energy.py` for independent closed
pilot/state review, `scripts/validation/collect_corrected_sampling_energy.py` for
the full 240-slot ledger, and `python -m redox.workflows.sampling_native --help`
for native preparation/dispatch. Only the three declared native pilots precede
the independently checked largest native G/resource receipt. Use
`scripts/validation/review_corrected_sampling_native.py` for native reconstruction
and `scripts/validation/observe_corrected_sampling.py` for
read-only process observations. Shared reservations, node preflight, aggregate
RAM/storage accounting, native thread caps and source-frozen namespace execution
remain in force. No old claim or output is cleared by age or an observation
timeout.

The saved-ionic-Hessian launcher correction exercised the electronic queue's
source-change guard. Its exact controller was confirmed closed before explicit
resume; existing workers and all results were preserved. The current source-bound
controller receipt is `electronic_queue/resume_launch_004.json` under the raw
energy root. This scheduling change does not alter the quantum method or QC.

The [six-receipt crosscheck](native_receipt_crosscheck_02.json) rechecks context,
source backing and all three component sums for every admitted native input. It
does not repeat physical quantum QC or complete an occupied ensemble. The saved
lead validation now has an actual exclusive H100 start. Our bulk coordinators
were explicitly stopped/resumed for the priority handoff while their workers
continued; the shared dispatcher honors a declared single-slot pending-lead
priority without changing methods or admitting an unqualified resource class.

The shared prelaunch accountant now separates recorded resource contention from
bounded environment/unclassified failures for the electronic, native and lead
controllers. It requires exact closed owners, absent scientific/guard outputs,
no launched-wrapper record, and the actual RuntimeError branch before every
Popen. Content-addressed frozen launcher versions preserve original traceback
line evidence after code updates. No claim or scientific output is modified.

The [audit](prelaunch_accounting_01.json) resolves five resource-refused inputs,
including three previously stranded at three receipts; their nonresource
failure counts are zero. Retry cooldowns remain, while resource refusals do
not consume the failure budget. Live/uncertain owners, started science and
unverified source branches remain held. The full 240-input objective and all
scientific QC gates are unchanged. [65 focused checks](shared_refusal_checks.xml)
cover ownership, source-version matching, resource accounting and terminal
cooldown behavior. All three coordinators were resumed after exact source-guard
closure, preserving their scientific workers. The [eight-receipt crosscheck](native_receipt_crosscheck_03.json)
retains every state's source/context and all three component sums.
