# Independent Nitro-pair review and gated NDI pilot

[Startup snapshot](summary.json) · [Tests](tests_dispatch_v3.xml) ·
[Pair waiter](../../../validation/uma_current_2026-10-24/H100_pair_cpu_preflight_2026-10-26/status.json) ·
[Live review observer](../../../validation/uma_current_2026-10-24/H100_pair_independent_observer_2026-10-26_v2/status.json) ·
[Live NDI pilot controller](../../../validation/ndi_external_uma_2026-10-26_v2/pilot_dispatch_v3/status.json)

An independent CPU observer now waits for the source-frozen fixed-device
FP32/Nitro comparison. When the pair finishes, it will recompute the energy and
force comparisons and check the declared inference modes, exact model identity,
same-device provenance and exclusive GPU residency. Failed, diagnostic or
incomplete evidence cannot pass. Even a passing result is limited to the frozen
fixed geometries; it does not certify optimizer behavior, DFT or experimental
accuracy.

The NDI neutral-state UMA resource pilot has a separate CPU coordinator. It
will not request a GPU until that independent pair review passes. It then
requires a free, at-least-32-GB GPU with the same device model and capability
as the parity comparison, fresh CPU-only preflight, atomic reservation, disjoint
CPU allocation and owned-process resource guards. Host and GPU memory are
sampled. An independent pilot review is required after completion; its result is
only a gas-phase warm start and resource measurement. The other two proposed
redox states are not automatically launched by this pilot.

At the [recorded snapshot](summary.json), both CPU controllers were live, the
pair was waiting for a suitable GPU, and the NDI pilot was waiting for parity;
neither had produced a new scientific result. Their live statuses can change.
Twenty-five focused tests cover parity rejection, device binding, preflight,
memory parsing, pilot review and isolated input routing. The earlier idle CPU
controller versions were stopped for source-provenance and preflight hardening;
no scientific worker was stopped.
