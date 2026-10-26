# Completed frozen local-motif interaction panel

[Reconstructed results](summary.json) ·
[Fragment construction](../local_motif_construction/summary.json) ·
[Earlier partial summary](../local_motif_interactions/summary.json) ·
[Construction figure](../../../figures/polymer_models/local_motif_construction.pdf)

All 36 frozen electronic states have independent reviews, and all 12 interaction
triplets now reconstruct at the original matched coordinates and settings. The
two close phthalimide control poses yield MeCN electronic interactions of
−24.71 to −24.81 and −26.20 to −26.24 kJ/mol across the three declared cap
lengths. The two separated poses are near zero: +0.047 to +0.049 and +0.060 to
+0.062 kJ/mol. These small positive values should not be overinterpreted.

Each value is `E(AB) − E(A) − E(B)` for two methyl-capped fragments in a frozen
parent pose, including cap and continuum response. The observed cap-length
spreads are under 0.10 kJ/mol for these four poses. This is not a conformer
ensemble, association free energy, decomposition of the full polymer cycle,
pure pi interaction, electronic coupling or an experimentally verified ranking.
Basis-superposition correction and broader basis sensitivity remain pending;
no gas-phase correction is transferred to the MeCN result.

The source-bound replay is
`calcs/validation/2026-10-26_local_motif/full_triplets_independent_replay_20261026.json`.
The existing collector used a fixed pending-status string even when all 12 sums
were present. The publisher checks all 36 source reviews, graph/coordinate
identities and inputs before issuing the explicit complete-panel status; the
original replay and prior partial summaries remain preserved. No new quantum
calculations were submitted for this publication step.
