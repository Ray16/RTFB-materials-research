# Guarded host switch-curvature qualification

The two later native switch-curvature tensors would occupy 43.37 GiB on the
115-atom site-3 geometry before other device storage. The preceding 32 GiB
pilot was stopped after that exact allocation lower bound was established; its
partial calculation remains unqualified.

The new profile partitions only the grid-point axis of the unchanged installed
analytic primitive. All atoms, switching weights, areas and point ownership
remain in each tile. Full tensors are held on the host and consumed in bounded
atom blocks. The native context's certified diagonal/off-diagonal curvature
bindings are refreshed when the transformed response functions execute, so a
namespace captured before native setup cannot silently restore dense primitives.

Both guarded full original-versus-modified spin controls pass:
[singlet](controls/1/status.json) and [doublet](controls/2/status.json).
Maximum Hessian differences are 7.86e-12 and 2.31e-12 atomic units. Energies,
coordinates and gradients also meet the unchanged tolerances. The guards retain
matching scientific PID/process-start identities, observe only the owned resident
on each assigned GPU, and report no interruptions. A solute-Hessian component
checkpoint is written before explicit solvent/CDS contributions; it is not a
qualified total Hessian, frequency set or G, and is not silently reused.

The native profile declares a 31 GiB measured-peak ceiling on a minimum 32 GiB
device, retaining one GiB physical headroom. The 28 GiB workspace value remains
an allocation estimate. This resource configuration does not change the native
force, symmetry, derivative, imaginary-mode, spin or partition criteria. Actual
large native resource/minimum/G qualification is still pending.

Sources: [grid-tiled host tensors](../../../scripts/pipeline/pcm_switch_curvature_blocks.py),
[scoped execution and checkpoints](../../../scripts/pipeline/smd_switch_curvature_workspace.py),
[full controls](../../../scripts/validation/qualify_switch_curvature_workspace.py),
[guarded same-process runner](../../../scripts/fleet/exclusive_python_worker.py)
and [native stages](../../../scripts/pipeline/switch_curvature_polymer_campaign.py).
The earlier response-only profile and its allocation-bound stop remain separate.
