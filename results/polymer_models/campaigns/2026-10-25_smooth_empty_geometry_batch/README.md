# Experimental-model empty-backbone geometry refinement

[Forty immutable plans](batch.json) cover twenty separated starts for each of
two declared backbone configurations. This is a new smooth-CDS model campaign;
no legacy energy or nuclear partition is transferred into it. The
[model controls](../../../validation/smooth_model_qualification_2026-10-25/README.md)
and [91-atom resource/geometry pilot](../../../reviews/2026-10-25/smooth_geometry_pilot/summary.json)
define its scope. A new final record uses `experimental_geometry_result.json`,
never the primary result schema.

[Dispatch](dispatch.json) is live orchestration state, not a completion
certificate. [Independent actual-start checks](../../../reviews/2026-10-25/smooth_batch_actual_start/summary.json)
observed seventeen resident scientific workers across lambda1, lambda2 and
lambda11, with disjoint CPU sets and native thread caps. Available suitable
GPUs are used through the shared reservation gate. A CPU coordinator's initial
SSH loader conflict was corrected before any scientific start; its
[recovery record](host_environment_recovery.json) preserves that distinction.

Every case has its own source/input-bound plan, gate owner, geometry progress
and durable raw data under
`calcs/polymer_models/2026-10-25_smooth_empty_geometry_batch/`.
Final graph, stereochemistry and multi-step force checks must pass. Geometry
completion supplies no active energy, Hessian, free energy, converged ensemble,
spacing rule or synthesis admission. The primary paired PMDI native calculation
remains a separate unchanged-model calculation.
