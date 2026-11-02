# Exact local viologen chloride-contact sensitivity

## Current checked thermochemistry

The initial N0-above stationary pose failed its native minimum test. The
[independent diagnosis](failed_minimum_diagnosis.json) reconstructs one unstable
internal mode at −50.52 cm⁻¹ from the preserved full Hessian, nine actual SCFs
and four sampled finite-difference checks. No frequency was made positive and
no G was admitted for that pose.

Both ±0.1-Å mode branches were optimized under stricter convergence and then
received entirely fresh electronic energies and native Hessians. Their
independent [minus](mode_minus_thermochemistry_review.json) and
[plus](mode_plus_thermochemistry_review.json) reviews pass, with lowest internal
frequencies 22.83 and 22.66 cm⁻¹. The independent partition reader validates the
same pinned archived certificate and source tree used by the actual worker;
it does not substitute the refactored reader's source hash.

The first finite-pose association cycles are now checked:
[minus branch](association_first_repaired_pose.json),
[plus branch](association_second_repaired_pose.json). At the 50-cm⁻¹ cutoff,
the 1-M values are +15.61 and +15.42 kJ/mol. The −7.92587-kJ/mol translation
conversion is applied once after the count-balanced source sum (association
changes solute count by −1). These are two repaired outcomes from one initial
pose, not independently weighted populations, an exhaustive conformer set,
an experimentally validated association constant or a grafting-yield prediction.
The five other checked contact starts still require full thermochemistry.

## Earlier preparation and geometry milestones

The [seven-start collection](geometry_collection_02.json) now contains six
independently checked stationary geometries, with the separated start still
running at that observation. Their primary r2SCAN-D4/SVP geometry-level energies
span 11.64 kJ/mol. The initial N0-above pilot is 4.72 kJ/mol above the lowest
checked geometry at that level. These are same-composition pose comparisons;
they supply neither association G nor independently weighted conformers.

The full primary electronic/native pilot has been prepared. Its first attempt
was [interrupted by foreign GPU residency](thermochemistry_interruption.json),
well below its host-memory bound, and supplies no qualified energy or G. The
fresh recovery preserves the original and transfers no qualification. Its first
placement was refused before science because declared host-RAM headroom was
unavailable; a fresh suitable reservation is required. All physical protocols,
resource limits and the exact input geometry remain unchanged.

The earlier [prepared association ledger](association_pending_02.json) keeps all
seven G cycles blank. Its maintained analysis verifies the exact primary parents and eventual
microstate component before matched arithmetic. The paired state is a single
statistical solute even though its graph contains disconnected chloride;
the association changes solute count by −1. A common base/HCl term is not added,
and the source translation-to-1-M correction follows the balanced sum. This
does not supply a complete salt/speciation model or experimental populations.

Shared actual-SCF checkpoint validation reproduces all nine existing molecular
native controls (81 frames) with unchanged scientific outputs; see the
[regression](primary_checkpoint_regression/summary.json). The
[focused checks](thermochemistry_checks_06.xml) pass 94 tests, including actual
CPU wavefunctions, scalar/frame/state corruption, standard-state/cluster-count
and incompatible-model rejection. Shared launch code records the exact budget
observer owner before observation starts.

The [experimental-source review](external_validation_evidence.json) distinguishes
qualitative molecular salt controls from quantitative association benchmarks.
No exact-graft MeCN association constant has been established. A separate
imidazolium/chloride source is retained as an unadmitted calibration candidate
until its full conditions, uncertainty and association definition are checked.

The [prepared manifest](prepared_manifest.json) contains seven geometry starts
for the exact oxidized local graft plus one chloride. All starts preserve the
covalent graft graph and contain a disconnected Cl⁻ component. The formal core
remains +2; the whole microstate has charge +1 and multiplicity 1. This is
chemically distinct from the isolated +1 doublet in the reduction ladder.

The parent is the checked `graft_ox2` structure in the corrected-primary CPU
example. Six deterministic chloride starts sample both nitrogen sites and
above/below/outward directions; a seventh begins 9 Å from the N–N midpoint.
Atomic separation and full hydrogen/heavy-atom connectivity are checked in the
actual XYZ files. These starting geometries are hypotheses, not measured
ion-pair populations or independent statistical weights.

The [first geometry review](pilot_geometry_review.json) now passes independent
execution closure, actual SCF/frame/occupied-spin, returned connectivity and a
fresh gradient check. Its maximum gradient is 2.8145e-5 Ha/Bohr and GPU peak is
10,610 MiB. This is a stationary geometry; a full Hessian is still needed to
qualify a minimum. The six remaining prepared starts are released within this
same measured geometry class, subject to live reservations and CPU/RAM/storage
capacity. All six now have [actual SCF starts verified](followup_start_verification.json),
with exclusive GPU residency, disjoint two-CPU bindings and all six native thread
caps set to two. The original [bounded dispatch](followup_dispatch.json) includes
one prelaunch GPU-contention refusal. Its exact stopped controller and absence
of science were checked before the [new placement](prelaunch_retry_dispatch.json);
the refused placement and receipts are preserved. No completed geometry or G
is inferred from these starts.

The [pilot plan](pilot_plan.json) retains the primary molecular SMD/MeCN model,
r2SCAN-D4/def2-SVP optimization, matched strict SCF/grid settings and actual
charge/spin. It does not substitute the separately scoped smooth-surface polymer
model. [Actual start](start_observation.json) confirms one live reserved GPU
worker on λ1, disjoint two-CPU affinity and entry into geometry evaluations.
Full SCFs, gradients and optimizer-returned coordinates are retained. The
initial start observation remains historical; it is not evidence of current
ownership. The later independent review releases the geometry class only.

No new association G, grafting yield, speciation population or material admission
is assigned. A later association cycle must compare this exact ion-associated
structure with the matched isolated dication and chloride under compatible
electronic/nuclear definitions. The change in solute count is −1, so standard-
state/concentration corrections cannot be copied from the count-balanced
separated-ion quaternization. Own Hessian/minimum, rotor, basis and final-state
checks remain required before a quantitative association conclusion.

Raw namespaces:
`/nfs/lambda_stor_01/homes/rzhu/redox_counterion_controls/2026-11-02_primary_chloride_contact_geometries_v2/`.
The first preparation is preserved separately; its parent manifest binding was
found to be hidden by the archived source overlay before launch, and is now
frozen inside each worker-visible source namespace. No source energies or G
values were transferred to the new microstate.

Use `python -m redox.workflows.counterion_geometry --help` for prepare/launch and
`scripts/validation/review_counterion_geometry.py` for independent closed geometry
review. The existing namespace preflight, atomic GPU gate, resource allocation,
full-state checkpointer and budget observer are reused. All cases remain scoped
as ion-convention controls for the established graft, rather than new synthesis
candidates.
