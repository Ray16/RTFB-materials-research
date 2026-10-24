# Exact paired-PMDI native thermochemistry

This fresh branch retains the checked 139-atom geometry, canonical optimized-state
gas/SMD energies and independently reconstructed fixed-density solvent-response
proof. [plan.json](plan.json) binds those exact parents and both full spin/Hessian
execution controls. The calculation uses SMD acetonitrile, the original quantum
levels and unchanged physical gates.

[launch.json](launch.json) and [startup_verification.json](startup_verification.json)
record the reserved V100, disjoint CPU allocation and actual scientific start.
The [native status](streamed_native/status.json) owns completion. SCF and
solute/explicit-PCM/CDS Hessian components are checkpointed on shared storage.
Missing native G remains blank until normal convergence, stationary geometry,
frequency, derivative, partition and source checks pass.

The execution profile streams solvent curvature, tiles weighted response and
suspends unused solvent matrices during pure-solute subroutines. The
[full controls](../../../../validation/solute_workspace_2026-10-24/summary.json)
verify exact matrix restoration and the complete original-versus-adapted Hessian.
The earlier full-native memory failure remains preserved in its own branch.

A [bounded CPU observer](../../../../reviews/2026-10-24/pmdi_pair_native_solute_followup/observer.json)
starts no quantum jobs. Only after successful native completion and verified
producer termination will it reconstruct the Hessian partition/geometry and
publish the checked component to the existing matched-cycle collector. Failure
preserves the current registry. No converged conformer ensemble or synthesis
admission is inferred from this one-minimum branch.
