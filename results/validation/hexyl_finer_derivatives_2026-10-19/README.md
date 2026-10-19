# Last-candidate finer derivative qualification

Both preserved hexyl-NDI radical geometries failed the coarse Hessian symmetry
check. The prior controller and numerical worker are terminal; no live job was
restarted or displaced. The new plan replays the independently qualified 1e-6-Bohr
CDS derivative with the same Hamiltonian, grids and physical acceptance thresholds.
Every geometry still needs native stationarity, minimum, symmetry and finite-gradient
checks, followed by active energies and independent admission review.

One suitable free H100 was available at submission. The launcher starts every free
suitable device up to the two pending inputs rather than requiring both devices
simultaneously. The [manifest](launch.json), [actual starts](start_verification.json)
and [resident-process check](residency_check.json) distinguish launched/running work
from completed output and the unscheduled second branch. The numerical step change
is qualification, not a repeated unchanged attempt. Maximum one actual attempt per
input; failures and old data remain preserved.

These large cases use the proven H100 memory class, disjoint two-CPU sets and native
thread caps. Derivative evidence stays frozen: active SPs are calculated in a separate
energy directory so no qualification input JSON is modified. Source optimizer
convergence is not inferred from old metadata; native physical checks and the
later independent review determine whether the preserved geometry is admissible.
Canonical data and primary predictions remain unchanged until that review succeeds.

## Terminal review

Both workers are now terminal and neither branch is admitted. The live PID and
GPU process probes confirmed exit. [Independent matrix/gradient checks](terminal_review.json)
reproduce the reasons for exclusion:

- Minus: the finer Hessian is symmetric, but its worst cavity-contact column
  differs from the finite-gradient derivative by 56.88%, above the 1% limit.
- Plus: the finer Hessian remains asymmetric by 2.3442e-5 Ha/Bohr², above the
  unchanged 1e-6 limit. No finite-difference/minimum success is inferred after
  this earlier gate failed.

All matrices, displaced gradients and failure records remain available. The
last candidate stays pending (24/25 checked potential pairs). Further progress
requires diagnosing the cavity-contact surface and finding a different
admissible minimum; repeating this same qualification is not scheduled.
