# Bounded PCM first-derivative prototype

Both full original-versus-modified spin controls passed, including complete
energy, optimized coordinates, gradients and Hessians. This prototype evaluates
only requested surface rectangles using the unchanged installed pairwise CUDA
primitive; the full off-diagonal derivative tensors are not allocated.

The [115-atom response pilot](../pmdi_site3_PCM_response_2026-10-23/README.md)
completed the blocked derivative stage but then exhausted memory constructing a
dense R operator. Its failure and raw baseline remain preserved. Small-control
success did not qualify the polymer resource footprint or native G.

The [checked successor](../pcm_first_response_workspace_v2_2026-10-23/README.md)
adds contracted R applications and a separately repeated pair of spin controls.
Producer sources, launch settings and small arrays are retained here; the large
failed response calculation is under `calcs/validation/pmdi_site3_PCM_response_2026-10-23/`.
