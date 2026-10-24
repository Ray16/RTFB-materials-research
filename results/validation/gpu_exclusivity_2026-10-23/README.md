# GPU residency guard and preserved device-sharing interruption

A post-launch audit found the two-graft redox optimizer sharing its H100 with a
foreign process. Only the exact owned worker was stopped, after preserving its
23 complete evaluated geometries. The peer process was left running. The initial
startup audit had missed the delayed CUDA attachment. [The incident record](../pcm_first_response_workspace_v2_2026-10-23/device_exclusivity_recovery.json)
records actual PIDs, GPU identity, the saved checkpoint and verified termination.
No molecular QC failure or accepted energy is inferred from this resource stop.

[watch_gpu_exclusivity.py](../../../scripts/fleet/watch_gpu_exclusivity.py) is a
CPU-only post-launch observer, complementing the authoritative `gpu_reserve`
allocation gate. It checks actual GPU UUID/PID residency every ten seconds.
Confirmed sharing stops only a matching same-user PID/process-start/campaign
command and its own descendants. It never signals peers, clears claims or GPU
locks, changes scientific settings or starts a calculation. Transient query or
ownership failures remain unresolved observations rather than restart triggers.

Finite observers are launched on lambda2 and lambda13, with separate validation/production observers on lambda1 through the original
candidate deadline. Each has one separately allocated CPU and native caps of one.
Their live owner/state files stay in the host subdirectories; launch receipts
record source hash, PID and deadline. Tests exercise actual process termination:
wrong identities are refused, the checkpoint is retained and a peer stays live.

The candidate observer also holds a phase with a stopped worker until an explicit
fresh recovery exists. The [new checkpoint continuation](../../polymer_models/campaigns/2026-10-23_crest_refinement/baseline__design_pmdi_ethyl__occ_2_3_exclusive_second_resume/README.md)
retains the old partial workspace, rechecks connectivity/stereochemistry and uses
its selected frame only as an unqualified warm start. Scientific outputs still
require the full optimizer, energy, force and thermochemistry checks.

`exclusive_python_worker.py` runs new short Python controls in the reserved process itself, retaining the gate parent and PID. Its five-second guard records actual residency during execution; guard files live outside the scientific tree being hash-frozen. The external observers continue to protect older production workers.

[The second interruption](second_residency_interruption.json) demonstrates the live guard catching another delayed foreign attachment. It stopped only the owned worker before an accepted result. Its partial attempt remains terminal; a fresh continuation uses the original evaluated frame and waits for an exclusive suitable device. No foreign process was signaled, and its later absence is not attributed to the guard.
