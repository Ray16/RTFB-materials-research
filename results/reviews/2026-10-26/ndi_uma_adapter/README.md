# NDI UMA adapter and Nitro comparison restart

[Checked startup snapshot](summary.json) · [Tests](tests_v2.xml) ·
[CPU-only preflight](preflight_snapshot.json) · [Ungated refusal](ungated_probe_v2.log) ·
[Prepared NDI campaign](../../../validation/ndi_external_uma_2026-10-26_v2/README.md)

The separate NDI campaign uses the existing UMA implementation through an
isolated manifest and seed directory. Tests with a mocked runner check exact
state routing, separate output paths, restoration after errors, and preservation
of a failed attempt. A real invocation without a GPU reservation was rejected
before creating a claim. These checks do not replace a real backend calculation.
The final adapter also refuses an interpreter mismatch instead of replacing a
process that could own a residency guard. The unstarted earlier draft is
[preserved](../../../validation/ndi_external_uma_2026-10-26/README.md).

The previous fixed-device FP32/Nitro comparison exhausted its wait window with
zero scientific starts. Its original records remain intact. A new CPU waiter
was started for the same frozen comparison, using a process-local substitution
of the CPU-only dependency preflight. At this snapshot its owner PID and start
time were live, the H100-host preflight passed without importing GPU modules,
and no validation worker had started. The waiter checks for an exclusive
reservation and a sampled quiet interval. Its live status is linked in
[summary.json](summary.json); a waiter is not a completed numerical comparison.

The fixed-geometry parity result, a real NDI UMA run, solution DFT and native
thermochemistry remain pending. The comparison does not by itself establish
optimizer behavior, throughput or experimental redox accuracy.
