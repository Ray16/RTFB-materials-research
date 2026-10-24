# Tiled weighted solvent response

Both full original-versus-adapted singlet/doublet controls pass, including
optimized coordinates, gas/MeCN energies, gradients and full analytic Hessians.
The maximum Hessian difference is 1.26e-11 atomic units, within the original
1e-6 limit. See [summary.json](summary.json) and the per-spin records under
`controls/`. Both workers ran on exclusively reserved V100s with disjoint
two-core allocations and matching native thread caps.

The new execution profile composes the existing streamed PCM curvature with
lazy diagonal weights and 256-row FP64 matrix products. It retains every atom
and surface point and the original quantum settings. This removes temporary
surface-square weighted matrices; the source library and earlier producers are
unchanged. Source and generated-function hashes accompany the evidence.

The [139-atom response test](../pmdi_pair_weighted_response_2026-10-24/plan.json)
uses the exact independently checked stationary geometry and canonical energy
parent. It passed and was [independently reconstructed](../pmdi_pair_weighted_response_2026-10-24/independent_review.json), with a maximum sampled relative error of 3.81e-6 and a 27,894 MiB GPU peak. Its [status](../pmdi_pair_weighted_response_2026-10-24/status.json) owns
completion. The fixed-density derivative check tests two response components at
0.0005 and 0.001 bohr; it does not qualify a full Hessian, minimum or free energy.
Raw outputs are in
`calcs/validation/pmdi_pair_weighted_response_2026-10-24/`.

The [previous response failure](../pmdi_pair_PCM_response_v2_2026-10-24/failure_review.json)
and its converged SCF checkpoint remain preserved. The two-graft native G and
matched PMDI occupancy cycle remain unresolved until the response and full
native stages pass their normal checks. No coupling, transfer rate, broader
conformer convergence or synthesis admission follows from these controls.

The first full native attempt failed in a separate pure-solute density-fitting allocation. The [successor profile and controls](../solute_workspace_2026-10-24/README.md) own the current native attempt, with SCF/Hessian checkpoints retained. The free energy and minimum remain pending.
