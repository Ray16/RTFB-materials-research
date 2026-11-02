# Complete corrected polymer-sampling cohort

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

The double-graft pilot stopped at its owned host-memory guard. Its main-process
RSS peak is 166.09 GiB but its summed process-tree peak is 332.03 GiB, above the
declared 240-GiB bound. Its last iteration checkpoint is retained unqualified.
Transient shared/forked memory is a hypothesis to test before recovery; the
record does not establish an out-of-memory event or permit ignoring the guard.
The single-graft pilot and the larger lead recovery continue separately. The
[later observation](progress_observation_02.json) records exact ownership and
eight completed response blocks for that lead at its observation time.

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
pilot/state review and `scripts/validation/observe_corrected_sampling.py` for
read-only process observations. Shared reservations, node preflight, aggregate
RAM/storage accounting, native thread caps and source-frozen namespace execution
remain in force. No old claim or output is cleared by age or an observation
timeout.
