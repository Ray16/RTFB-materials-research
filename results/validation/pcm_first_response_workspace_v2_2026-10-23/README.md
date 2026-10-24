# Bounded PCM first derivatives and contracted Hessian R

This execution profile retains the installed pairwise CUDA primitive, surface
points, float64 precision and direct LU solver. Requested derivative rectangles
are assembled from unions of point blocks. Three Hessian-response functions use
the existing solvent R operator rather than materializing a dense R matrix.
Installed library files are unchanged.

Full original-versus-modified [singlet](controls/1/status.json) and
[doublet](controls/2/status.json) controls pass. Maximum Hessian differences are
7.82e-12 and 1.86e-12 atomic units; energies, coordinates and gradients also meet
the original tolerances. Execution receipts confirm that every transformed
response path was exercised. CPU tests independently check dense primitive
values, overlapping/empty slices, transpose signs, translation invariance,
asymmetric R batch orientation and fail-closed template checks.

The [exact 115-atom response pilot](../pmdi_site3_PCM_response_v2_2026-10-23/README.md)
passes independent fixed-density matrix finite differences, with a maximum
relative error of 4.15e-6 and measured GPU peak 31,224 MiB. It qualifies only the
stated response/resource test, not a full native Hessian, minimum, G or polymer
ensemble. The 28 GiB workspace estimate is an allocation heuristic; the isolated
response test's declared physical ceiling is 32 GiB. The subsequent native pilot
retains its own separate resource ceiling and all physical checks.

A [fresh native campaign](../../polymer_models/campaigns/2026-10-23_crest_refinement/baseline__design_pmdi_ethyl__occ_3_native_first_response/README.md)
requires the exact geometry's response proof and imported energy identity. It
still runs the normal native resource/physical qualification; preparation does
not establish that the complete workflow fits a V100.

Sources: [first-derivative blocks](../../../scripts/pipeline/pcm_first_derivative_blocks.py),
[derivative adapter](../../../scripts/pipeline/smd_first_response_workspace.py),
[contracted-R adapter](../../../scripts/pipeline/smd_first_response_workspace_v2.py),
[full controls](../../../scripts/validation/qualify_first_response_workspace_v2.py)
and [native stages](../../../scripts/pipeline/first_response_polymer_campaign.py).
The first prototype and its failed large test remain separate.

The two-graft optimizer subsequently had a device-exclusivity interruption.
[The incident record](device_exclusivity_recovery.json) preserves its 23 evaluated
checkpoints and verifies that only the owned worker was stopped. The
[residency guard](../gpu_exclusivity_2026-10-23/README.md) now watches subsequent
production jobs; no stopped checkpoint is promoted to an electronic result.

The first full V100 native attempt was subsequently stopped after inspecting two later switch-curvature allocations: their 43.37 GiB lower bound alone exceeds the device. [The bound record](native_allocation_bound.json) and old partial calculation are retained. The [guarded host-switch successor](../pcm_switch_curvature_guarded_2026-10-23/README.md) addresses that additional stage; response-only success remains separate from full-native qualification.
