# Strict-matmul Nitro diagnostic

This separate first-case test retains the FP16-off control and sets the installed executor's internal `_USE_TF32` policy false in the reserved process. It records the original/effective policy, every FP16 switch and source hashes. The internal setting is restored on exit. This is an implementation diagnostic, not an official default profile or production qualification.

The [launch record](launch.json) declares a finite reserved H100 test using the original model checkpoint, input frames, reference engine and repeatability thresholds. Live status is in `dispatch/status.json`, with raw calls under `calcs/uma_nitro_validation/2026-10-27_ieee_intermediates/`. No scientific start or outcome is inferred from preparing the waiter. The result requires independent review and cannot qualify the full ten-case panel by itself.

The test has now ended. [Independent reconstruction](independent_review.json)
confirms that reference force repeatability passes (1.68×10⁻⁶ eV/Å), while the
internal strict-matmul Nitro control fails (1.96×10⁻⁴ eV/Å) against the unchanged
1×10⁻⁴ eV/Å threshold. Both energy-repeat checks pass. Each engine retained 29
calls, including eight calls on the identical frame, on the same H100 and model
checkpoint. The maximum later Nitro repeat difference is 2.39×10⁻⁴ eV/Å.
This does not establish a unique implementation cause. Neither this internal
profile nor the default Nitro profile is qualified for the current production
panel; changing the threshold is not inferred from the result.
