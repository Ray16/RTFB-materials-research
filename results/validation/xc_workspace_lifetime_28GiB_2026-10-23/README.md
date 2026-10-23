# 28 GiB XC execution controls for the V100 native pilot

Both complete original-versus-modified small implementation controls passed:
[singlet](controls/1/status.json) and [doublet](controls/2/status.json).
The process-local adapter changes allocation order and completed AO-batch lifetime
with unchanged equations, grids, precision and scientific gates. Its separate
[139-atom gradient qualification](../xc_workspace_lifetime_2026-10-23/README.md)
is a geometry-level resource/derivative check, not a native-Hessian certificate.

The checked site-3 electronic state was handed to a
[separate native-thermal campaign](../../polymer_models/campaigns/2026-10-23_crest_refinement/baseline__design_pmdi_ethyl__occ_3_native_xc_surface_budget/README.md).
The V100 worker has an actual reserved start; no accepted native G exists yet.
The preflight verifies dependencies, load, memory and disjoint two-core affinity
with native thread pools capped at two.

A geometry-specific host forecast uses the measured maximum cavity size of
24,615 points, plus 1,024 safety points, in the original bounded-PCM memory
formula. The 128 GiB cache limit and 128 GiB extra headroom are unchanged. The
resulting 216.16 GiB forecast fits the recorded available RAM and 75% admission
fraction. The earlier blanket 256 GiB admission did not fit that fraction.
This changes the declared resource estimate, not the chemical model or QC.

Earlier unlaunched preparations and their producer bodies are preserved under
`calcs/validation/xc_workspace_lifetime_28GiB_2026-10-23/prelaunch_source/`.
Raw controls and pilot checkpoints stay under the corresponding `calcs/` paths;
[launch.json](launch.json) and the native campaign plan retain their identities.
