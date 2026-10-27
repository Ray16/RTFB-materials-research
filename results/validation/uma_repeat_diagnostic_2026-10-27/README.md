# Retained Nitro fixed-geometry repeatability failure

The [independent reconstruction](independent_review.json) reads every saved call
from a same-H100 comparison of the exact first failed PMDI neutral case. The
standard FP32 reference and Nitro use the same checkpoint and four frozen input
frames, with the original repeatability tolerances. Each mode has 29 calls,
including eight evaluations of the identical first geometry.

The reference's initial-to-final force error is 1.64984e-6 eV/Angstrom; Nitro's
is 7.57206e-4, exceeding the unchanged 1e-4 limit. Both energy repeat errors
are below the declared 4.2e-4 eV bound. Among later repeated first-frame calls,
the largest Nitro force difference is 7.83648e-4 eV/Angstrom, so the retained
variation is not confined to the first call. The implementation cause is not
established by this observation.

The full-panel comparison and NDI pilot remain unqualified. A separate
[precision control](../uma_repeat_fp32_control_2026-10-27/README.md) tests the
installed executor's explicit FP16-intermediate switches without changing the
production profile or tolerance. Saved reference/Nitro calls and outcomes are
under `calcs/uma_nitro_validation/2026-10-27_first_case_repeat/attempt_001/`;
launch, preflight, gate and residency evidence remain in `dispatch/`.
