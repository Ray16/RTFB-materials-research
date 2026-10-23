# XC temporary lifetime: implementation and fixed-checkpoint checks

The previous 139-atom two-graft optimizer stopped after 23 durable energy/gradient
evaluations when the GPU XC gradient routine had no available workspace. Its
entire trajectory is retained. A completed density loop still held its last AO
batch while entering the next workspace stage. The controlled process-local
adapter releases that completed batch and synchronizes before the next stage;
the equations, numerical grids, precision and physical acceptance gates remain
unchanged. No installed library was edited.

Both original-versus-modified small spin controls passed under the 64 GiB
execution profile: [singlet](controls/1/status.json) and
[doublet](controls/2/status.json). The launcher observation timed out while the
actual calculations continued and completed; this was not a scientific failure
and no duplicate control was submitted.

The preserved evaluated checkpoint 23 was checked separately on a V100 using a
28 GiB execution budget. The [point report](checkpoint_139_v100/gradient_qualification.json)
and [independent reconstruction](independent_checkpoint_review.json) verify
sampled total-energy derivatives for two components at two steps (eight displaced
SCFs plus the baseline). Maximum derivative error is 1.27e-7 Ha/Bohr versus the
1e-6 limit; the maximum gradient difference from the stored parent is 8.08e-9
Ha/Bohr. Observed peak device memory is 27,808 MiB. The initial checkpoint graph
and stereochemistry were also rechecked.

This qualifies a geometry-level gradient test on 32 GB hardware. It does not
qualify the full hybrid-energy/Hessian workflow, prove an optimized minimum or
establish conformer convergence. A [fresh full continuation](../../polymer_models/campaigns/2026-10-23_crest_refinement/baseline__design_pmdi_ethyl__occ_2_3_xc_checkpoint_resume/README.md)
uses the saved coordinates only as an unqualified warm start, with all ordinary
optimization and physical checks still required.

Producers: [XC adapter](../../../scripts/pipeline/xc_workspace_execution.py),
[spin controls](../../../scripts/validation/qualify_xc_workspace.py),
[checkpoint preparation](../../../scripts/validation/qualify_xc_checkpoint.py)
and [full-resume preparation](../../../scripts/pipeline/prepare_xc_checkpoint_resume.py).
Plans and receipts bind inputs, controls and producers by hash. Raw diagnostics
are retained in `calcs/validation/xc_workspace_lifetime_2026-10-23/`.
