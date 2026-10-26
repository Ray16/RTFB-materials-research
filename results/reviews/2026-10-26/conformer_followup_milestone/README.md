# Broader conformer follow-up

This is an execution and software-validation checkpoint, not a new polymer
free-energy result. [Actual worker inspection](actual_workers.json) records
project GPU residency, exact process start identities, disjoint CPU allocations
and native thread caps on three nodes. Scientific completion requires the
independent stage reviews below.

The two recovered occupied-graph pilots now support the final
[38 geometry inputs](../../../polymer_models/campaigns/2026-10-26_smooth_remaining_occupied_expansion/manifest.json).
Together with the earlier 152 occupied inputs, ten checked occupied references
and forty empty-backbone inputs, the declared refinement set contains twenty
starts for each of twelve graph/configuration sets. Starts are not minima or
Boltzmann weights. The [coverage audit](coverage.json) checks all 240 exact input
hashes, distinct raw namespaces and ten occupied-class resource reviews.
[Independent review of the additional geometries](../remaining_occupied_geometry/status.json)
keeps failed or unfinished calculations outside accepted results.

Separate [single-graft](../../../polymer_models/campaigns/2026-10-26_smooth_occupied_SP_single_capacity/batch.json)
and [double-graft](../../../polymer_models/campaigns/2026-10-26_smooth_occupied_SP_pair_capacity/batch.json)
hybrid-energy pilots measure their own resource requirements. Their prospective
256-GiB host ceiling is a bound, not measured resource evidence. They use the
same independently checked exact geometry for gas and MeCN calculations; no
energy is transferred from the empty chain. [Their independent reviewer](../occupied_electronic_pilots/status.json)
runs separately.

The [empty-chain native pilot](../../../polymer_models/campaigns/2026-10-26_smooth_native_pilot_checked/batch.json)
retains each full displaced gradient and SCF record, the Hessian and the nuclear
partition. [Independent reconstruction](../smooth_native_pilot/status.json)
checks the sampled Hessian columns, internal modes, detected symmetry number,
25/50/100 cm^-1 thermochemistry and actual resource/ownership records. The
[source-bound follow-up request](../../../polymer_models/campaigns/2026-10-26_smooth_native_followup/request.json)
can declare forty exact-parent native states only after this pilot passes.
The existing reference may be reused only for the identical geometry; the other
39 states require new Hessians and native free energies. Electronic RAM-recovery
options remain attached to the same exact geometry.

The new native driver and dependency/resource checks pass 38 focused tests;
the [test report](native_capacity_tests.xml) records the cases. Negative cases
include NaN/infinite or inconsistent free energies, a mismatched model, state,
basis, runtime or parent geometry, incomplete electronic review, resource-bound
violations and inappropriate promotion to primary results. This is separate from
the previously checked full software suite.

Raw calculations are under `calcs/polymer_models/` at the exact namespaces in
each plan. Launch/CPU-affinity recovery records are under
`calcs/validation/2026-10-26_followup_launch/`. The first unstarted CPU launches
inherited the launcher's narrow affinity and were rejected before any scientific
attempt; separate attempts corrected the startup mask so every child allocates
its own CPU slot. Original logs and ownership records remain intact.

The occupied energy placement also exposed a RAM-preflight mismatch before any
scientific attempt: the dependency probe required 64 GiB, although these pilots
require their 256-GiB ceiling plus 32 GiB headroom. The 128-GiB node therefore
rejected placement before an attempt directory or worker existed. All four
intended hosts were inspected; only the declared read-only observer matched.
[The recovery](occupied_capacity_recovery.json) preserves the original dispatch
and intent records, applies the actual bound during CPU preflight, and treats a
late capacity shortage as contention before spawning. The native continuation
uses its own 80+32 GiB requirement. Only its CPU controller was replaced; the
live native quantum worker and every physical input remain unchanged.

The same placement adapter now serves the empty-chain electronic stages.
[Five unstarted cases](empty_energy_capacity_recovery.json) were inspected on all
four submission hosts before their intent records were archived and reopened.
The two 80-GiB cases' earlier RPC cause remains unattributed; no attempt directory,
scientific stage, GPU guard or matching producer existed. The three 128-GiB RAM
experiments also require more memory than the 128-GiB node can supply with
headroom. Thirty-five checked electronic records and the three original
80-GiB scientific failures remain intact. Only the verified CPU controllers were
replaced; each physical worker and exact-coordinate energy input is unchanged.

All these stages use the separately identified smooth-surface MeCN model. They
cannot replace primary-model PMDI or phthalimide records. Broader ensemble
convergence, motif-specific close-contact energetics and electronic coupling
remain unresolved. The checked manuscript continues to report finite sampled
cycles with these limitations.
