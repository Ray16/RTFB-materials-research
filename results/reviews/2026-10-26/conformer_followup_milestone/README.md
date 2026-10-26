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

Separate [single-graft](../../../polymer_models/campaigns/2026-10-26_smooth_occupied_SP_single/batch.json)
and [double-graft](../../../polymer_models/campaigns/2026-10-26_smooth_occupied_SP_pair/batch.json)
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

The new native driver and dependency/resource checks pass 34 focused tests;
the [test report](native_dependency_tests.xml) records the cases. Negative cases
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

All these stages use the separately identified smooth-surface MeCN model. They
cannot replace primary-model PMDI or phthalimide records. Broader ensemble
convergence, motif-specific close-contact energetics and electronic coupling
remain unresolved. The checked manuscript continues to report finite sampled
cycles with these limitations.
