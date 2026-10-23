# Two adjacent PMDI grafts: full workspace continuation

The exact 139-atom resource point now passes with a peak of 39,206 MiB and both
hybrid SCFs converged. Its original OOM outputs remain unchanged. This separate
full campaign retains the same graph, initial geometry, numerical model and QC;
its original pilot evidence is explicitly reused, rather than reported as a new
calculation. Fresh force storage avoids the earlier pilot/optimization collision.

[plan.json](plan.json) records sources and inputs. The standard
[worker record](energy/workers/lambda13_gpu0_queue.json) identifies the original
optimizer, now stopped after an XC-gradient workspace allocation failure. The [execution guide](../../../../../docs/POLYMER_WORKSPACE_EXECUTION.md)
describes resource-profile reuse and native-stage handoffs. The candidate
observer now follows the separately prepared XC checkpoint successor.

Raw calculations and durable geometry evaluations:
`calcs/polymer_models/2026-10-23_crest_refinement/baseline__design_pmdi_ethyl__occ_2_3_workspace_64GB/`.
Neither the starting point nor intermediate optimizer evaluations supply a
stationary-state free energy, conformer convergence or synthesis admission.

The evaluated trajectory is preserved. The [XC checkpoint successor](../baseline__design_pmdi_ethyl__occ_2_3_xc_checkpoint_resume/README.md) uses checked starting coordinates and unchanged final acceptance requirements; no intermediate value from this failed optimizer is promoted to a minimum or G.
