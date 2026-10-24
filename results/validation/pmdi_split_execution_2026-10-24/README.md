# Separate polymer geometry and active-energy execution

This execution branch preserves the established SMD(acetonitrile) geometry and
active electronic-energy protocols while storing their stages independently.
It addresses the resource needs of the exact 139-atom, five-repeat, two-graft
PMDI graph. It does not change a functional, basis, solvent, precision or
physical-admission threshold to fit the device.

The [frozen pilot plan](plan.json) contains two separate reference points.
The [GPU geometry point](geometry/status.json) passes on a 32 GiB V100:
27,898 MiB peak memory, zero electronic-energy difference and
3.64e-8 Ha/bohr difference in the maximum gradient versus the original preserved
checkpoint. This is a fixed-point resource/reference check, not an optimized
minimum or a native Hessian qualification.

The original [CPU active-energy point](electronic/status.json) used the unchanged
ωB97M-V/def2-TZVPD gas/SMD calculations at the original large-device reference
geometry. It exceeded the 300 GiB RSS bound during CPU PCM integral allocation
and stopped with no accepted energy. Its complete resource trace and failure
remain preserved. [Separately checked PCM tile-budget controls](../cpu_pcm_memory_2026-10-24/README.md)
now support a fresh large reference-point run. Its energies must agree within
1e-5 eV, and convergence, host-memory and scratch limits must all pass. Running
or partial SCF evidence supplies no accepted electronic energy. The new CPU
pilot must finish before the optimized geometry can use that execution route.

The [separate full geometry stage](../../polymer_models/campaigns/2026-10-24_split_execution/baseline__design_pmdi_ethyl__occ_2_3_geometry/plan.json)
has an actual reserved V100 optimizer start. It begins at the preserved geometry
checkpoint and writes every new optimizer evaluation and SCF recovery record to
durable storage. Final graph/stereochemistry, total-energy derivative and
stationarity checks remain required. It deliberately produces neither an active
hybrid energy nor native G.

Launch accounting reserves the running site-2 Hessian's **remaining** host-memory
bound; its existing RSS is already excluded from available RAM. The geometry
worker has a separate 192 GiB RSS ceiling and 32 GiB free-RAM floor. Its finite
host and GPU-residency guards signal only their owned scientific process.
[Two subprocess checks](host_guard_tests.xml) verify normal same-PID execution
and an interrupted worker with preserved resource evidence.

Scripts: [resource pilots](../../../scripts/validation/qualify_split_polymer.py),
[geometry stage](../../../scripts/pipeline/optimize_polymer_geometry.py), and
[host guard](../../../scripts/fleet/bounded_host_python_worker.py).
The [CPU energy handoff](../../../scripts/pipeline/finish_polymer_electronic.py)
requires both completed qualifications before preparing a new calculation. It
uses the canonical active-SP writer, retains both geometry and energy parents,
and produces an ordinary electronic completion receipt. Its result can be
adopted by the native stage without another hybrid calculation. [Nine focused
checks](handoff_tests.xml) include mismatched parents and a forbidden-recompute
test. These software checks do not qualify an unfinished CPU or geometry pilot.
Raw geometry states live under
`calcs/polymer_models/2026-10-24_split_execution/`; pilot outputs and launch,
input and source hashes are retained beside this guide. Native minimum/G,
matched-cycle reconstruction, sampling convergence and candidate decisions
remain separate requirements.
