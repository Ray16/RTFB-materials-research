# Nitro intermediate-precision diagnostic

This separate first-case comparison disables ten installed executor FP16
intermediate switches before importing the executor, then records their actual
runtime values and source hash. It uses the original frozen first-case inputs,
checkpoint, reference engine and repeatability thresholds. No production profile
or full-panel qualification is changed.

The [launch record](launch.json) identifies a finite reserved H100 diagnostic.
Live preflight, placement and outcome are under `dispatch/`; raw calls are under
`calcs/uma_nitro_validation/2026-10-27_fp32_intermediates/`. A spawned waiter is
not a scientific start. The result requires inspection of the recorded calls,
effective flags and retained tolerances before any interpretation or further
comparison. This diagnostic alone cannot qualify all motifs or geometries.
