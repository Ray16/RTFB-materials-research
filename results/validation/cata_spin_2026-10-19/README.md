# Cata-NDI dianion alternate-spin comparison

The source state's active SP became available, so the required DFT alternative
was launched under its own system reservation. The [plan](plan.json) freezes its
exact chosen geometry, charge, multiplicity, active protocol and energy-record
hash before any GPU is claimed. The [launch](launch.json) and [follow-up observation](../cata_solution_2026-10-19/followup_start_verification.json)
record actual resource ownership and progress separately from comparison results.

The maintained `redox.validation.spins` workflow performs a full-response
alternative geometry and matching active SP. Near electronic gaps require
native spin thermals; near free-energy gaps retain ensemble-admission and kinetic
gates. A lower alternative can be promoted only after full physical checks and
exclusive molecule writer ownership. Existing native state writers retain their
locks. Failed helpers and gate-limited comparisons stay visible. No unverified
spin energy enters a primary candidate potential.

The calculation is complete: the triplet is higher by 0.158664 eV electronically
and 0.141495 eV in matched native free energy. The reviewed two-spin model gives
a 0.4041% triplet population and a -0.000104032-eV equilibrium correction. An
[independent policy review](ensemble_review.json) passes 38 checks.
[Selection](selection_receipt.json) is waiting for the remaining molecule writer;
the prepared model is not yet a primary correction. Component spin degeneracy is
counted once. No experimental solution-rate validation is inferred.
