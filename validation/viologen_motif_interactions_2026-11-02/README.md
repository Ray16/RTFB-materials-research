# Frozen viologen motif interactions

The four checked charged-polymer parent frames now have complete A/B/AB
electronic controls at three cap lengths. All twelve triplets and 72 actual
SCF checkpoints pass independent frame, occupied-orbital spin, source, model,
resource and arithmetic review. The
[comparison](motif_definition_comparison.json) and
[complete ledger](motif_energy_collection.json) keep these electronic quantities
separate from association G, coupling and whole-chain nonadditivity.

| Parent | Whole-chain ΔΔE | Frozen motif E interaction, MeCN | Frozen motif E interaction, gas |
| --- | ---: | ---: | ---: |
| All-same, adjacent | −3.46 | +9.03 | +335.72 |
| All-same, one spacer | +18.68 | +7.62 | +283.91 |
| Alternating, adjacent | +1.20 | +9.63 | +357.45 |
| Alternating, one spacer | −13.48 | +7.79 | +289.93 |

All energies are kJ/mol at the same r2SCAN-D4/def2-SVP level and corrected
implementation. The MeCN model is the separately declared experimental smooth
surface used by the ionic polymer calculations. Each value represents one
retained parent frame. The 1.07/1.09/1.11-Å cap variation changes MeCN motif
interactions by at most 0.00265 kJ/mol; this is a deterministic cap-length
sensitivity, not an uncertainty estimate or chemical-cap/basis convergence.

The cut retains each benzyl carbon, replacing its aryl bond by hydrogen to form
a methyl-capped viologen dication. Every retained coordinate stays fixed. The
28-atom fragments each have charge +2 and singlet multiplicity; their frozen
56-atom pair has charge +4 and singlet multiplicity. A CH2-constrained attachment
query prevents the symmetric viologen core from selecting its other N-methyl
carbon. The archived constructor and independent pose/cap auditor are reused.

The two definitions have opposite signs in two of four parent frames. A negative
whole-chain occupancy cycle therefore does not demonstrate attractive local
motifs under this cap model. Its difference from the frozen interaction is not
assigned uniquely to backbone, linker, caps, deformation or solvent. No binding
free energy, counterpoise correction, electronic coupling, experimental
population, graft loading or synthesis priority is inferred.

Raw construction and calculations remain on durable shared storage:

- `/nfs/lambda_stor_01/homes/rzhu/redox_motif_interactions/2026-11-02_viologen_frozen_methyl_caps/`
- `/nfs/lambda_stor_01/homes/rzhu/redox_motif_interactions/2026-11-02_viologen_methyl_energy_pilot_v2/`
- `/nfs/lambda_stor_01/homes/rzhu/redox_motif_interactions/2026-11-02_viologen_methyl_energy_followups/`

The earlier pilot preparation, which caught a missing source-path binding before
submission, is retained separately. The versioned preparation, calculation and
review entry points are `viologen_motif_fragments`, `motif_interaction_energy`
and `review_motif_interaction_energy.py`. This folder contains self-contained
SCF scalar/state receipts and the exact frozen XYZ inputs for manuscript replay;
that replay does not repeat the original physical QC or fresh QM.

## Checked presentation and manuscript milestone

The immutable `2026-11-02_three_ionic_cycles` manuscript epoch now integrates
this electronic counterexample and three independently checked finite ionic G
cycles. The article/SI build at 20/33 pages. The slide deck has 16 main and six
supporting slides; all pages were visually inspected, with detailed checks of
new slides 13 and 20. Source hashes, decimal cycle/motif replay, generated TeX
and compiler checks pass. All 228 indexed local result paths and the result
README links resolve (`result_links_check.json`). The self-contained local
reviewer export also compiles and matches both checked documents. No Overleaf
remote is connected.

The fourth ionic G cycle remains blank after a checkpoint-file I/O failure;
its Hessian and completed evaluations are preserved. Occupied ensembles, broader
counterion/speciation controls, material eligibility, public deposition and
author metadata remain open. This is a checked presentation/draft milestone,
not a submission-ready manuscript.
