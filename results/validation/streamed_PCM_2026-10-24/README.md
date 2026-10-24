# Streamed solvent-Hessian execution on V100

The [full singlet and doublet replays](replay_qualification.json) pass against
the unchanged dense calculation on reserved 32 GiB V100 GPUs. The maximum
full-Hessian differences are 1.38e-11 and 1.07e-12 Ha/bohr², respectively.
Energy, geometry and gradient comparisons also pass the original bounds.
Both controls exercise every streamed operator. These small controls establish
numerical equivalence for their inputs, not a large-polymer free energy.

The [implementation](../../../scripts/pipeline/streamed_pcm_execution.py)
contracts analytic switch-curvature tiles directly and uses bounded GPU
off-diagonal curvature blocks. It retains FP64, all surface points, the
functional, basis, solvent and physical checks. The GPU block cache is capped
at 128 MiB and each primitive union at 512 MiB. Installed libraries are unchanged.

The [large site-2 recovery](../../polymer_models/campaigns/2026-10-24_streamed_native/baseline__design_pmdi_ethyl__occ_2/plan.json)
has an [actual scientific start](../../polymer_models/campaigns/2026-10-24_streamed_native/baseline__design_pmdi_ethyl__occ_2/startup_verification.json)
on lambda1 GPU 0, with two disjoint CPUs, an 80 GiB RSS bound and a 32 GiB
available-RAM floor. This is a running resource and physical qualification;
its native G and matched occupancy cycle remain pending. The established
stationary electronic parent is reused without recomputing its hybrid energy.

The [native driver](../../../scripts/pipeline/run_streamed_native.py) uses the
existing native-DFT thermochemistry producer and finer-CDS derivative checks.
It preserves separate solute, explicit-PCM and CDS Hessian components before
combination, plus partial switch contractions. A [same-process wrapper](../../../scripts/fleet/checkpoint_scf_worker.py)
adds durable SCF density/orbital files, coordinates, settings and cycle energies.
Partial components and SCF files are unqualified until the complete physical
and resource checks pass. Final independent reconstruction remains required.
Raw native records live under
`calcs/polymer_models/2026-10-24_streamed_native/`.

The [previous worker review](../pmdi_site2_runtime_2026-10-24/summary.json)
records an owned stop after prolonged CPU solvent-curvature work. Its
interruption traceback localizes the slow path to CPU `get_d2D_d2S`; it does
not demonstrate a geometry failure or GPU hardware fault. Its partial solute
Hessian and original records remain preserved and are not substituted into G.

[Twenty-five focused checks](combined_tests.xml) pass, including tensor-axis
and transpose identities, signed-weight contractions, allocation refusal,
SCF callback preservation and UMA state/parity contracts. Software checks and
small spin replays do not establish conformer convergence or candidate eligibility.
The [full regression](regression_tests.xml) passes 1,044 tests. Its two
resource/dependency skips each pass separately with [two CPUs](two_CPU_test.xml)
and the [figure environment](figure_dependency_test.xml), respectively.
The [milestone check record](milestone_checks.json) also verifies control and
input hashes and all 128 local guide links.

For UMA inference, use the separate [V100 runtime qualification](../uma_v100_2026-10-24/README.md).
V100 uses standard inference; actual Nitro requires newer GPU architectures.
