# Expanded finite chloride-contact sensitivity

The [current seven-outcome collection](current_03.json) incorporates the
[N0-outer plus branch](N0_outer_plus_native_review.json), whose
[matched association cycle](association_N0_outer_plus.json) is +5.93 kJ/mol
at 1 M and 50 cm-1. Seven outcomes from five original starts form **five
numerical duplicate groups** at every tested threshold. The N0-above and
N0-outer signed pairs each share one group, with proper-rotation heavy-atom
RMSDs of 0.00377 and 0.01472 Angstrom. All 21 outcome/cutoff values are positive;
the primary range remains **+5.33 to +22.36 kJ/mol**. Individual G values remain
visible, and neither pair supplies hit-count entropy or population weights.

The article, SI Table S11 and supporting slide 21 include all seven outcomes
under the new immutable evidence epoch. The maintained collector is
`python -m redox.validation.contact_cohort --help`, with the explicit
[roster](seven_outcome_roster_02.json). It binds the original geometry campaign,
exact native reviews and matched ledgers. Shared outcome identities and
geometry-derived grouping replace assumptions that outcome count determines
the number of starts or duplicate groups. Historical four/six-outcome snapshots
remain supported. [Fifty checks](seven_publication_checks_03.xml) cover both old
cohorts, the new cohort, numerical/state/context changes, wrong branch ledgers
and reordered rosters. The [committed source](source_export.json) also passes
the [full isolated manuscript validation build](committed_source_validation.json).
No population, synthesis or material inference is added.

## Earlier six-outcome milestone

Two additional contact components pass independent full thermochemistry review:
[N0-below, minus repair](N0_below_minus_review.json) and
[N0-outer, minus repair](N0_outer_minus_review.json). Their lowest internal
frequencies are 24.45 and 20.96 cm-1, and sampled-column relative errors are
6.15e-5 and 1.41e-4. Actual electronic and nine native SCF checkpoints, spin,
geometry, Hessian, modes, nuclear terms and closed execution are checked.

Their matched primary-model 1-M association G values at the 50-cm-1 rotor cutoff
are +6.23 and +5.33 kJ/mol. The [six-outcome collection](current.json), retaining
all four earlier outcomes, spans **+5.33 to +22.36 kJ/mol**. All eighteen
outcome/cutoff combinations remain positive. The [independent Decimal replay](independent_arithmetic.json)
reconstructs each paired-minus-graft-minus-chloride sum and applies the
two-to-one, 1-atm-to-1-M correction once (-7.92587 kJ/mol); maximum discrepancy
is below 4e-10 kJ/mol.

Proper-rotation, symmetry-aware heavy-atom comparisons give five numerical
duplicate groups at all three declared RMSD thresholds, with the same
0.2-kJ/mol electronic tolerance. Only the two earlier N0-above repairs merge.
Six search outcomes are not six statistical weights. These finite primary
MeCN results do not establish complete ion populations, iodide exchange,
experimental association accuracy or grafting yield. The covalent
quaternization cycle is unchanged.

The earlier four-outcome publication snapshot remains immutable. The new
contacts preserve the positive finite-association interpretation but lower
the sampled minimum by about 10.09 kJ/mol. All six outcomes are now integrated
into the [3 November manuscript evidence epoch](../../manuscript/evidence/2026-11-03_six_chloride_contacts/manifest.json),
the article, SI Table S11 and supporting slide 21. The source-bound
[publication summary](publication_summary.json) and
[38 publication checks](publication_checks.xml) document this update; the
earlier editorial snapshot is preserved.

The [N0-outer plus geometry](N0_outer_plus_review.json) also passes independent
review. Its own thermochemistry is now prepared and queued under the existing
32-GB-class envelope. The [new manifest](followup_manifest.json) retains the
three earlier attempts, including the unqualified N0-below plus minimum
failure. The [controller handoff](controller_resume.json) proves the previous
owner closed; no started scientific attempt was retried. All original failures
and snapshots remain available. `collect.py` records the source-bound union
and duplicate analysis used here.
