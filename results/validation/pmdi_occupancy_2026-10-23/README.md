# Exact PMDI candidate: polymer feasibility and capacity scenarios

This campaign follows `design_pmdi_ethyl`, exactly matching enumerated design
`search_pmdi_77642a20bc89`. Its capped-molecule redox results are preserved. The
[precursor review](../../reviews/2026-10-23/pmdi_precursor_access/summary.json)
documents access to the configured mono-N-ethyl/NH precursor. The complete route
still has status **hold**; these feasibility models do not admit synthesis or a
new full-redox campaign.

Three five-repeat graphs have grafts at site 2, site 3, or both adjacent sites.
They contain 115, 115 and 139 atoms. Unoccupied reactive sites retain benzyl
chloride; positions 1 and 5 are styrene padding. The open backbone has ten carbons,
nine backbone C–C bonds and two ends in every actual XYZ and drawing. An identical
empty backbone is identified as a shared component; occupied phthalimide energies
are never transferred to PMDI.

- [Frozen source-search plan](plan.json), [candidate snapshot](candidate_snapshot.json) and [route snapshot](route_snapshot.json).
- [Exact identity and geometry input review](input_validation.json).
- [Actual source-search starts](progress_01.json): three 14-core CPU searches in ALPB acetonitrile.
- [Conditional capacity scenarios](capacity_scenarios.json): formula/graph-checked upper bounds.
- [Baseline DFT declaration](baseline_plan.json), [current components/cycle](baseline_current.json), and [finite observer](baseline_coordinator.json).
- [Construction figure](../../figures/polymer_models/pmdi_occupancy_construction.pdf) and [capacity-bound figure](../../figures/polymer_models/pmdi_capacity_bounds.pdf).

Validated MMFF starts also enter a separate, limited DFT baseline while broader
CREST searches run. Each baseline case has a real energy resource pilot, full
geometry/energy checks and a native thermal pilot. A checked single minimum does
not establish a converged conformer ensemble. Case directories share the common
resource root so launchers account for existing RAM commitments. Missing
electronic energies, G values and the conditional cycle remain blank.

The 139-atom resource pilot exhausted GPU memory in SMD matrix assembly. The two smaller pilots completed their fixed-geometry evaluations but exceeded the unchanged 80 GiB safety budget. A [control-qualified execution branch](../smd_build_memory_2026-10-23/README.md)
reduced temporary-buffer lifetimes, but that larger case later exhausted memory in a solvent-response product;
[the recovery registry](baseline_recoveries.json) preserves the original failure
and links the separate resource requalification. Numerical equivalence on the
controls does not by itself establish larger-case memory suitability.

At three handles per five repeat units, the assumed two-electron upper bounds
for one, two and three grafts are 61.2, 98.9 and 124.5 mAh/g. They retain scaffold
and residual-chloride mass, exclude unlisted crosslinker/electrolyte mass, and
assume full utilization. The one-electron alternative is also reported. These
are deterministic composition scenarios, not inferred loading probabilities or
measured accessible capacity. A target above a bound rejects that scenario;
a target below it does not establish material feasibility. No target was chosen
from the computed ranks.

Raw searches: `calcs/polymer_models/pmdi_occupancy_2026-10-23/`.
Raw DFT/native baseline: the same root's `baseline/` subdirectory.
The [machine entry](../../polymer_models/candidate_occupancy/current.json) links
these outputs into normal candidate publication. Backbone configuration, finite
fragment mass and sampling scope remain explicit.


The [bounded-workspace successor](../smd_bounded_memory_2026-10-23/README.md)
passes full original-versus-modified spin controls. Its same-model three-case
recovery registry keeps the original 80 GiB budget and all physical gates.
Available larger-case fixed-geometry energies/gradients must also agree with
the preserved references. New large pilots await a free suitable H100.
Raw recovery locations are given by each current branch's `plan.json` under
`calcs/polymer_models/2026-10-23_crest_refinement/`; the earlier generic raw-path
statement does not describe these recovery branches.


The later [32 GB gradient recovery](../smd_gradient_workspace_32GB_2026-10-23/README.md) passes the exact site-2 point and sampled energy derivatives. A [full geometry/energy branch](../../polymer_models/campaigns/2026-10-23_crest_refinement/baseline__design_pmdi_ethyl__occ_2_workspace_32GB/README.md) has started on lambda2. The candidate-only observer linked above is historical during this controlled handoff; its exact PID was verified and stopped without stopping scientific workers. The wider-search observer and all fifteen source searches remain live at the recorded check. The other occupied baseline components and thermal stages remain pending. This continuation does not change the current all-axis decision.
