# Empty-backbone sampling in both configurations

[Summary](summary.json) · [Complete analysis](analysis.json) ·
[Figure](figure.pdf) · [Tests](tests.xml) · [Current native census](native_census.json)

Both twenty-input empty-backbone sets have checked native thermochemistry in
the separate smooth-surface acetonitrile model. Full-atom connectivity, hydrogen
attachments and configured stereochemistry are checked before heavy-atom
alignment. Symmetry mappings preserve chirality and use proper rotations;
reflections and silently truncated symmetry enumeration are rejected.

At RMSD cutoffs of 0.05, 0.10 and 0.20 Å and a 0.20 kJ/mol electronic-energy
duplicate tolerance, the all-same set contains 20 distinct structures and the
alternating set contains 19. Alternating inputs 16 and 17 belong to the same
duplicate group. Each well contributes once; search-hit counts are not weights.
The groups defined on the full twenty-input set are restricted to the first ten
for a consistent nested comparison.

| Configuration | Minimum E change | Minimum G change | Ensemble G change |
| --- | ---: | ---: | ---: |
| All same | 0.00 | 0.00 | −1.03 |
| Alternating | −0.20 | 0.00 | −1.23 |

Changes are 20 minus 10 inputs, in kJ/mol, at 298.15 K, a 50 cm⁻¹ rotor cutoff
and a 0.10 Å duplicate cutoff. All minimum-E, minimum-G and ensemble-G changes
remain within the declared 2 kJ/mol criterion across all nine rotor/duplicate
settings. Ensemble changes span −1.07 to −0.85 kJ/mol for all same and −1.28 to
−1.13 kJ/mol for alternating. These are procedural sensitivities, not confidence
intervals or a bound on undiscovered minima.

The configurations are analyzed separately and are not Boltzmann pooled into
an assumed atactic-polymer population. Both empty components meet this finite
input-batch diagnostic. This does not close the occupied ensemble cycles or
prove global conformer convergence. Source searches and selected inputs retain
their original completion/sampling limits.

The analysis binds forty exact native/energy parent pairs and their scored
geometries. Three additional native reconstructions, including the last two
empty states and a separate site-4 resource point, are preserved in
[native_additions.json](native_additions.json). Together with the preceding
replays, the working census is 75 checked native states; the additional site-4
point is not included in the empty-backbone partitions.

The maintained analysis is
`scripts/validation/analyze_native_conformer_batches.py`, with reusable geometry
and partition functions in `src/redox/properties/conformer_geometry.py` and
`conformer_ensemble.py`. A separate standard-library manuscript checker
reconstructs duplicate grouping, finite partitions, settings comparisons and
the numerical LaTeX macros from frozen records. Thirty focused checks include
atom-order invariance, rejection of reflected stereocenters, duplicate-hit
entropy, changed energies, incorrect populations and unsupported pooling.

Raw replay inputs are listed in [inputs.json](inputs.json). Quantum calculations
remain under `/nfs/lambda_stor_01/homes/rzhu/redox_polymer_continuations/` and the
unchanged archive. No new quantum job runs in this analysis. Primary-model
PMDI/phthalimide cycle components are not replaced with these alternate-model
values.
