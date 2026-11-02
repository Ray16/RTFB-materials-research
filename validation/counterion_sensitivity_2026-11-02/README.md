# Exact local viologen chloride-contact sensitivity

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
