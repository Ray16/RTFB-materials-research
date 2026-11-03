# Expanded finite chloride-contact sensitivity

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
the sampled minimum by about 10.09 kJ/mol. They therefore warrant inclusion
in the next scientific evidence update before submission; the finalized
editorial PDFs have not been silently replaced.

The [N0-outer plus geometry](N0_outer_plus_review.json) also passes independent
review. Its own thermochemistry is now prepared and queued under the existing
32-GB-class envelope. The [new manifest](followup_manifest.json) retains the
three earlier attempts, including the unqualified N0-below plus minimum
failure. The [controller handoff](controller_resume.json) proves the previous
owner closed; no started scientific attempt was retried. All original failures
and snapshots remain available. `collect.py` records the source-bound union
and duplicate analysis used here.
