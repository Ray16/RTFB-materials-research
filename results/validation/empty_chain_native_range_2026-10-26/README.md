# Empty-chain conformer free-energy range

[Results](summary.json) · [Frozen reviewed inputs](manifest.json) · [Tests](tests.xml)

Four independently reviewed native free energies are currently available for
the same five-repeat, empty benzyl-chloride backbone and one specified backbone
configuration, in a separately identified smooth-surface acetonitrile model.
At the 50 cm⁻¹ rotor cutoff their range is **7.19 kJ/mol**; the ranges at 25
and 100 cm⁻¹ are **7.96** and **7.42 kJ/mol**, respectively. The lowest of these
four is the second declared input. These values illustrate the sensitivity of
an empty-chain reference to the sampled geometry.

This is a range among four checked records, not an uncertainty interval or a
converged ensemble. The calculation does not establish that all four optimized
geometries are distinct basins. It cannot quantify the error in a grafting cycle,
compare tacticities, characterize an interaction, or determine experimental
loading. The alternate solvent model must not be mixed with the primary-model
cycle. The original reviews contain the physical-QC and raw-record provenance;
this diagnostic binds their hashes but does not repeat their underlying quantum
calculations or all raw-record checks.

Reproduce from the repository root with
`python scripts/validation/summarize_checked_empty_conformers.py --manifest results/validation/empty_chain_native_range_2026-10-26/manifest.json --root . --out /path/to/new-summary.json`.
The input list is frozen; later native records need a separate manifest and
review, not silent substitution into this snapshot.
