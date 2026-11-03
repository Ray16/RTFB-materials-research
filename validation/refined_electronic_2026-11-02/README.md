# First electronic scoring at accepted phenothiazine refinements

All three refined states now pass independent electronic review:
[neutral minus](phenothiazine_neu_mode_minus_review.json),
[neutral plus](phenothiazine_neu_mode_plus_review.json) and
[oxidized doublet](phenothiazine_ox_stationarity_review.json). Ten actual SCF
checkpoint frames support both implementation variants. Corrected-minus-installed
differences are at most 9.10e-13 Ha; each selected occupied wavefunction passes
spin and orthonormality review. No old-geometry energy reference was inserted.
Their [own native preparations](native_prepared_manifest.json) and
[three archive preflights](native_preflight.json) now pass; the bounded native
coordinator waits for suitable free 32-GB-class capacity. G and potential remain
blank until each own Hessian, minimum, derivative and nuclear review passes.
The [first own-native start](native_actual_starts_01.json) is now independently
observed on a V100 with exact worker ownership, disjoint CPUs, required thread
caps, exclusive residency and its live budget observer. Its native center SCF
has begun; the other two inputs retain their queued slots. This start does not
qualify the Hessian, minimum or G.

The [current ledger](current.json) retains three exact inputs: both neutral
singlet branches and the oxidized +1 doublet. All three geometries independently
pass their strict final-gradient, graph/hydrogen, actual SCF-frame and occupied
spin checks. They describe N-(4-methylbenzyl)phenothiazine, a local graft-site
model rather than a full polymer. The same reviewed exact-route gate remains
required. This does not admit the held F/CN proposals.

Each new geometry needs its own primary omegaB97M-V/def2-TZVPD/SMD MeCN energy.
The maintained installed/corrected electronic worker and its original solver,
multiple-guess selection, checkpoint, convergence and resource checks are reused.
The [prepared inputs](prepared_manifest_v2.json) have 39 atoms and 976 high-level
basis functions. A new explicit profile leaves `reference_energy_Ha` null and
both stored-energy comparisons blank. Only the installed-versus-corrected
difference can be tested here; old-geometry energies are never substituted for
an absent reference. The independent reviewer reconstructs the actual selected
occupied-wavefunction spin as well as frames, settings and scalar sums. The
solver's four-decimal spin diagnostic is compared at its recorded precision;
full-precision reconstruction and the existing spin limit are retained.

The [archive preflight](archive_preflight_v2.json) checks all three source-visible
namespaces without starting CUDA. The copied admission profile is also executed
inside the actual archive mount. The [first actual start](actual_scientific_starts_01.json)
verifies the neutral-minus worker on a 32-GB V100, exact host/PID/start, two
disjoint scientific CPUs, all six required runtime caps of two, exclusive GPU
residency and a separate live budget observer. Actual SCF settings have been
written. This is startup evidence, not a completed or qualified electronic
energy. Pending E, G and potential fields remain blank.

The shared bounded coordinator also queues three previously checked chloride
geometries for their own electronic/native components. Their [preparations](../counterion_T4_geometry_2026-11-02/followups/thermochemistry_prepared_manifest.json)
retain the validated full 32-GB-class Hessian envelope; the T4 geometry
qualification cannot transfer. Both bulk sampling coordinators were resumed
after their exact old owners exited on the shared launcher's source guard.
The [current owner observation](controllers_actual_owners.json) confirms all
four active admission processes. Their scientific workers retain immutable
calculation sources and were not signaled. All eligible free GPUs can be used;
occupied devices and unqualified T4 Hessian capacity remain unavailable.

An [initial CLI adapter failure](initial_adapter_failure.json) returned before
the existing budget observer: the underlying launcher returns a process/record
pair. Every original wrapper/observer owner was proven closed, with no
scientific program, guard, SCF or energy output. Original namespaces and claims
remain preserved. The corrected adapter keeps the budget observer alive until
its owned child closes. Fresh v2 namespaces supersede those unused inputs;
no started scientific calculation was retried. The [112 focused tests](final_checks.xml)
cover this lifecycle regression, exact state/resource admission, historical
reference refusal, the two maintained dispatcher adapters, contact physical
context and existing queue controls. A [fresh legacy replay](legacy_oxidized_review_01.json)
also reconstructs the original oxidized state's actual checkpoints and all
three historical comparisons unchanged. No new quantum calculation is inferred
from that replay.

The native continuation uses one shared constructor for original and refined
primary components. It preserves the new XYZ bytes and exact electronic parent,
retains the source optimizer record, and changes only the local result's ID/state
labels needed by the native engine. The numerical worker is unchanged. Both
neutral branches and the oxidized state require 21 actual center/displaced SCFs,
the standard two columns plus sulfur x/y/z at both steps, full Hessian/mode and
own nuclear/spin reconstruction. [128 final checks](native_final_checks.xml)
cover both native profiles and all preceding launch/scope contracts.

An additional [spin diagnostic](oxidized_spin_diagnostic_01.json) identifies three
converged oxidized guesses per variant. One lies only 4.4e-10 Ha above the true
minimum but rounds its spin contamination differently. The verifier now selects
the exact lowest-energy checkpoint before comparing the rounded producer
diagnostic. Its actual S-squared is 0.76774763 and the reported 0.7677 is correct.
The physical spin threshold is unchanged, the raw calculations are preserved,
and a regression rejects substitution of a nearby guess. As-used reviewer code
is retained byte-identically in [the source index](as_used_review_source_index.json).

Use `python -m redox.workflows.route_refined_electronic --help` for prepare,
launch and bounded queue entry points. Review a completed placement with
`python scripts/validation/review_screening_energy.py review --help`.
Use `python -m redox.workflows.counterion_thermochemistry queue --help` for the
separate full-contact queue. The first-scoring campaign supplies no native G:
each accepted geometry still needs its own full solution Hessian, minimum and
sampled-derivative review, sulfur-centered columns, nuclear partition and
context-matched Fc reference before an oxidation potential can be assembled.
Both neutral branches remain explicit until those checks support comparison.

The first native-minus component subsequently failed its unchanged symmetry
gate. The [independent saved-matrix diagnosis](neutral_minus_Hessian_symmetry_failure_02.json)
rechecks the closed physical execution, actual center SCF, finite 117-by-117
matrix and 6.543e-6 Ha/Bohr-squared asymmetry (limit 1e-6). The maximum pair
is C11-x/H32-x. Sampled derivatives, minimum and G are unqualified. The
[versioned read-only diagnosis](diagnose_saved_hessian.py) preserves the matrix
and settings; this scientific attempt was not retried or symmetrized to pass.
The partner branch continues in its own namespace.

Three closed incoming-chloride GPU refusals were initially counted as generic
prelaunch errors because the shared reader recognized only the screening
launcher. The [reconstruction](contact_resource_refusals_reconstructed.json)
now verifies the exact molecular-delegate traceback before Popen: all three
are resource waits, with zero scientific attempts and zero failure-budget
consumption. Shared accounting snapshots both actual launcher sources and
keeps live/unknown owners, started wrappers and guards held. [181 tests](resource_delegate_checks.xml)
cover both delegates and the existing scientific gates. The
[exact controller handoff](resource_accounting_closed.json) and
[resumed owner observation](resource_resumed_actual_owners.json) record four
coordinators on disjoint CPUs; in-progress scientific workers were not signaled.
These latest records supersede earlier owner/start observations above.
