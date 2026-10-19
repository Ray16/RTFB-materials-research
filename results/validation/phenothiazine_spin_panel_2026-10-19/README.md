# Alternate-spin checks for the experimental parent panel

The [plan](plan.json) binds seven ready primary states from three exact parents.
The [launch](launch.json) and [actual start verification](start_verification.json)
record three reserved V100 workers on lambda1, with disjoint two-CPU sets and
all native thread pools capped at two. One worker owns each ready molecule and
checks its states sequentially, respecting the existing molecule spin lock.
Pending primary states hold no GPU in this spin campaign.

Each source-addressed check optimizes the alternate singlet/triplet or
doublet/quartet surface and compares active-protocol SMD energies. Near or lower
alternatives require their own qualified native minima and corrected thermal
partitions. Lower-spin installation waits for the primary molecule writer to
finish. Relevant mixtures remain blocked until their explicit equilibrium policy
is qualified; successful electronic comparison alone does not prove experimental
accuracy, global conformer populations or kinetics.

`parents/<id>/status.json` records progress and the case reports. Detailed immutable
case evidence is under the paths returned by `redox.validation.spins`; alternate
geometries are under `calcs/spincheck/<id>/`. There is one real attempt per selected
case, with preserved failures and the existing three-failure circuit breaker.
New primary states need a fresh ready-input batch after their current molecule
worker releases its spin lock; they are listed in `pending_primary_states`.

The helpers use the corrected optimizer-convergence and complete-cache transaction
paths, qualified by 466 software tests. Keep plan-bound numerical producer files
unchanged while these workers are active. A new code epoch requires separate
frozen sources and a new run plan; do not relabel existing evidence.
