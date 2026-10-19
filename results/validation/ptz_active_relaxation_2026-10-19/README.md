# Active-level phenothiazine geometry sensitivity

The checked free-parent DMeOEPT second oxidation differs from experiment by
+0.547 V. Its reference-free charging gap is 0.9531 V calculated versus 0.599 V
measured. Prior SCF seed checks found no lower sampled electronic solution at the
frozen geometry. Active-level forces at the composite geometries are about
0.015 Ha/Bohr and pass finite-energy controls. A local lower-level-Hessian estimate
suggests only -0.0084 eV of charging-gap change, but is not an actual optimization.

This [three-state plan](plan.json) tests that estimate directly. The neutral,
radical cation and dication relax independently at the same active wB97M-V/TZVPD,
SMD/MeCN electronic level, with full grid/NLC response. The SCF settings are the
previously checked tighter gradient-control settings. Optimizer and final SCF
convergence, declared graph, cold final energy, stationary gradient and a worst
component finite-energy check are all required. Checkpoints and failures remain
preserved; each state has one actual attempt.

The [node preflight](lambda13_preflight.json), [launch](launch.json) and
[actual optimizer starts](start_verification.json) record three suitable reserved
H100s, disjoint CPU pairs and capped pools. This is a private sensitivity model:
the geometry treatment changes while the electronic level is fixed. No canonical
geometry, energy, thermal record, benchmark residual or candidate export changes.
A stationary geometry alone does not prove a Hessian minimum or a consistent new
free energy. No fitted offset or correction is permitted.

After all three states finish, independently reconstruct
`E(neutral) + E(dication) - 2 E(radical cation)` and compare its change against the
original electronic gap. Ferrocene cancels. Keep the old thermal gap separate;
this job does not recompute thermochemistry or establish a replacement potential.
The [requirement audit](../../reviews/2026-10-19/requirement_completion_audit/audit.json)
records why the original goal remains active despite the passed software and
arithmetic checks. The published candidate workbook/slides remain authoritative.
