# Bounded vinyl-monomer solution pilot

This run covers three exact vinyl-bearing monomers: NDI ethyl, PMDI ethyl and
ethyl viologen, with three declared charge states each. It was requested after
the broad autonomous goal was paused; that broader goal remains paused.

Six reserved H100 GPUs on lambda13 started the finite nine-state queue. The
[launch receipt](launch_lambda13.json), [node preflight](lambda13_preflight.json)
and [actual start verification](start_verification.json) distinguish submissions
from resident predictors. The initial verification found six actual DFT workers,
zero completed states and zero failed states. These counts are an as-of record.

Each state screens three force-field conformers with xTB/ALPB warm starts,
r2SCAN-D4/SVP(D)/SMD optimization and wB97M-V/TZVPD/SMD active energies, followed
by the existing derivative-checked native Hessian/internal-mode thermochemistry.
Optimization stopping tolerances are tightened; the energy surface and physical
admission thresholds are unchanged. The [plan](plan.json) freezes inputs and
numerical producers. One real attempt is permitted per state; failed outputs are
preserved for review. Claims record host, PID and process-start identity.

Outputs live under `workspace/calcs/dft/<vinyl_id>/<state>/`. Task reports and
stage logs live under `tasks/<vinyl_id>__<state>/`. Completion receipts are in
`done/`; failures are in `failures/`. Exact private structures are never mixed
with capped-model records, and the shared scientific registry is unchanged.

A completed state receipt verifies that declared spin surface. Independent
review, ground-spin selection and complete-couple checks remain required before
publishing potentials or adding numerical results to the slide deck. The
current figure collection retains the pending labels until that admission.

Run the read-only observer on the worker host with the absolute redox interpreter:
`python observe.py`. It verifies process ownership, GPU residency, disjoint CPU
sets, thread caps, claims and actual stage-log progress. Inspect existing owners
and preserved attempts before any restart; this pilot does not auto-retry failures.
