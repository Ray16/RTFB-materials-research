# Empty-backbone 10-to-20 conformer comparison

[Comparison](comparison.json) · [Ten-input manifest](manifest_10.json) ·
[Twenty-input manifest](manifest_20.json) · [Finite twenty-state result](finite_20.json) ·
[Independent arithmetic](independent_arithmetic.json) · [Tests](tests.xml)

The first twenty search-ranked empty-backbone conformers in the `all_same`
configuration have independently checked native free energies in the separately
identified smooth-surface acetonitrile model. The twenty structures remain distinct
under all three heavy-atom RMSD duplicate thresholds (0.05, 0.10 and 0.20 Å),
combined with the 0.20 kJ/mol electronic-energy criterion. Each minimum contributes
once; search-hit frequency is not a population weight.

Expanding the first ten to the first twenty changes neither the lowest electronic
energy nor the lowest individual free energy. The finite-set ensemble G changes
by **−1.03 kJ/mol** at 298.15 K and the 50 cm⁻¹ rotor cutoff. Changes at 25 and
100 cm⁻¹ are −0.85 and −1.07 kJ/mol. All nine combinations of duplicate threshold
and rotor cutoff meet the previously declared 2 kJ/mol batch-stability threshold.
The twenty-state configurational lowering relative to the lowest individual G
is 4.52 kJ/mol at 50 cm⁻¹.

This establishes the declared batch diagnostic for **one empty graph/configuration**.
It does not prove global sampling convergence, complete the occupied-state cycles,
or validate experimental conformer populations. The alternate-model values cannot
be mixed with primary-model grafting or PMDI cycles. Approximate qRRHO and explicit
conformer entropy may overlap for soft modes. The other backbone configuration and
the occupied graphs remain separate requirements.

The producer, `scripts/validation/compare_empty_conformer_batches.py`, reconstructs
both stored finite analyses, verifies their source-bound native/electronic records,
checks exact graph/configuration/method identity and requires ranks 1–10 and 1–20.
The geometry comparison uses `analyze_empty_conformer_set.py`; partition arithmetic
uses `src/redox/properties/conformer_ensemble.py`. Raw paths and hashes are in the
finite-set inputs and bound native reviews under `calcs/polymer_models/`. A separate
stdlib log-sum reconstruction agrees. Twenty-seven focused tests pass, including
rejection of changed stored partitions and postselected or incompatible batches.
No new quantum calculations were launched by this analysis.
