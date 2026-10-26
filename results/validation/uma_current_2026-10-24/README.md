# Current default UMA entry point

`redox uma` and the normal `redox launch` now prefer Nitro on supported GPUs.
V100 selects the retained standard calculator explicitly. The
[installation check](entrypoint_installation.json) confirms the installed public
command and unchanged versions of all 239 distributions. See the
[execution guide](../../../docs/UMA_BACKENDS.md) and [policy](../../../config/uma_execution.json).

The state-aware calculator invalidates ASE results when charge or spin changes
at identical coordinates, creates a new fixed-state predictor and propagates
explicit cache resets. Its tests cover spin scans, false Nitro fallback,
graph/reference-bound warm starts, checkpoint status and shard/CPU ownership.
Current-model warm starts have their own cache and identity; no model-ID
override or legacy-energy equivalence is assumed.

The [actual V100 comparison](legacy_comparison.json) passes all ten states and
forty geometries: maximum differences are 2.35e-6 eV and 1.74e-6 eV/Å.
[Cached points](legacy_points.json), [launch](legacy_launch.json) and
[entry-point tests](entrypoint_tests.xml) retain the evidence. This validates
the fallback integration, not Nitro execution on V100 or a Nitro speedup.

The previous [fixed-device official FP32/Nitro waiter](H100_pair_quiet/status.json)
ended its finite window with no validation start. A [new CPU-preflight waiter](../../reviews/2026-10-26/ndi_uma_adapter/README.md)
has been started for the same frozen comparison and is tracked separately.
Reservation contention consumes no scientific attempt. The comparison uses the
original frozen ten-state/four-frame bounds and checks actual inference modes.
A started failure remains for review rather than being repeatedly resampled.
Raw paired results will be written under
`calcs/uma_nitro_validation/2026-10-24_current_model_pair_quiet/attempt_*/`.

All solution-energy, stationarity, native-Hessian and candidate-admission checks
remain necessary. Execution preference is not a numerical migration certificate.
