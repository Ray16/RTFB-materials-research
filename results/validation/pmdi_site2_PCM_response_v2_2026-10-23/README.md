# Exact site-2 response/resource qualification

[The declared test](plan.json) uses the separately qualified saved site-2 geometry.
[Completed matrices](status.json) and [independent reconstruction](independent_review.json)
pass all four fixed-density comparisons at 0.0005/0.001 bohr. Maximum relative
error is 4.978e-6; observed GPU peak is 26,552 MiB on the reserved V100. A CPU-only
residency observer monitors the registered process without introducing a second
CUDA context.

The analytic response has shape 115 × 3 × 1150 × 231. It holds the baseline AO
density and MO coefficients fixed; displaced evaluations are electrostatic
potential calculations, not new SCF or native free energies. Raw matrices,
geometry, SCF checkpoint and source hashes remain under
`calcs/validation/pmdi_site2_PCM_response_v2_2026-10-23/`.

This independently qualifies the response component for this exact geometry.
The later full native pilot still requires its own device/host suitability,
Hessian, minimum, spin and nuclear-partition checks. No polymer G or complete
occupancy cycle is inferred from the response test.
