# Solvent derivative validation and polymer continuation

This is a numerical-method checkpoint, not a polymer-completion certificate.
Use [independent workflow validation](../../polymer_models/workflow_validation.json)
and [live completion status](../../polymer_models/completion/status.json) for later progress.
The solvent remains SMD acetonitrile and all physical admission thresholds remain unchanged.

| Question | Evidence | Result |
| --- | --- | --- |
| Does refining the electronic grid fix the open-trimer force? | [Paired same-coordinate calculations](../polymer_gradient_common_axes_2026-10-21/plan.json) | No; both grids show about 1.81e-4 Ha/Bohr error at atom 7/Y |
| Which force term is responsible? | [CDS energy derivatives](../cds_energy_gradient_2026-10-21/summary.json) | The legacy CDS term reproduces the discrepancy; conservative derivatives match full-energy differences |
| Are conservative derivatives step converged? | [Smaller-step cross-check](../cds_energy_gradient_finer_2026-10-21/summary.json) | Open trimer passes; the old stacked checkpoint needs further review |
| Are spacing geometry forces affected? | [Frozen-source six-case audit](../spacing_CDS_forces_frozen_2026-10-21/summary.json) | All six selected geometries agree with CDS energy derivatives |
| Where does the spacing Hessian asymmetry arise? | [Component diagnostic](../polymer_hessian_components_2026-10-21/status.json) | Default CDS step 1e-4 Bohr; refined 1e-5 step reduces total asymmetry from 2.06e-6 to 2.05e-8 |
| Does additive reconstruction reproduce a complete Hessian? | [Two full-recomputation replays](../hessian_reconstruction_replay_2026-10-21/summary.json) | Singlet/radical maximum matrix errors 1.34e-11/2.83e-9; thermal differences below 1.2e-8 eV |

The [spacing continuation plan](../../polymer_models/campaigns/2026-10-21_spacing_cds_recovery/plan.json)
pins source matrices, additive CDS components, both complete reference replays and
all six geometry-force audits. Its largest-model pilot reuses the unchanged
electronic/PCM curvature, replaces only the refined CDS component, and repeats
native SCF, force, finite-gradient, spin, minimum and qRRHO checks. Additional
two-step finite-gradient checks target the worst CDS columns. The pilot has now passed: maximum asymmetry 2.05e-8, lowest vibrational
frequency +7.23 cm−1 and standard finite-gradient errors below 7.2e-5 relative.
The remaining four missing Hessians have actual running processes on separate
H100s; [resource verification](../../polymer_models/campaigns/2026-10-21_spacing_cds_recovery/parallel_start_verification.json)
checks their GPU residency, disjoint CPU sets and thread caps. Together with the
previous empty-backbone record, two native free energies qualify.

The first [open-trimer restart](../../polymer_models/campaigns/2026-10-21_pmdi_open_conservative/plan.json)
passed its initial conservative force check, then stopped at the next geometry
when the two-step derivatives differed by 1.17e-5 Ha/Bohr. That failed attempt
remains intact. Separate [open](../../polymer_models/campaigns/2026-10-21_pmdi_open_adaptive/plan.json)
and [stacked](../../polymer_models/campaigns/2026-10-21_pmdi_stacked_adaptive/plan.json)
continuations refine unresolved CDS derivative components without changing the
1e-6 Ha/Bohr acceptance tolerance. Each derivative input is saved before evaluation.
Both deterministically perturbed stacked seeds pass the initial derivative audit;
the declared positive perturbation is the optimization seed, without energy selection.
These two continuations are prepared and waiting for suitable unoccupied GPUs. The final geometry must
also pass the original optimizer, topology/stereochemistry, active-SP and
independent gradient checks, plus full-energy derivative checks on its two
largest force components at two steps. No Hessian minimum is inferred from an
electronic receipt.

All historical failures, checkpoints and matrices remain at their original
shared-storage paths. The [checkpoint inventory](summary.json) hashes diagnostic
evidence and records the new raw roots. The [regression results](software_tests_final.xml)
contain 751 passed tests and one figure-dependency skip; that rendering test is
checked separately in the figure environment.

The current single-conformer spacing analysis includes all twelve electronic
conformers. Its selected electronic nonadditivities are −10.98 and −14.50 kJ/mol
for adjacent and one-spacer scenarios. This change from the earlier incomplete
panel demonstrates sensitivity to the reference conformer; it is not evidence
of a universal spacing trend. Conditional free-energy cycles remain blank until
the exact selected geometries' thermal records qualify. These short, neutral
method controls do not establish bead association or electron-transfer rates.

The [shared-filesystem restart experiment](../shared_restart_2026-10-21/summary.json)
uses 32 actual design recipes, deduplicated to 17 exact molecular tasks. A lambda1
worker exits abruptly after four checked tasks. Recovery verifies its death on
lambda1; lambda2 executes the remaining 13 without repeating those four. All 17
resumed descriptor records match uninterrupted execution and serial definitions.
This qualifies that descriptor stage, not quantum restart equivalence.
