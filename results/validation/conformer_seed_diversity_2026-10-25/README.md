# Geometry diversity of conformer inputs

The [inherited-input audit](summary.json) finds that twenty saved starts need
not represent twenty separated geometries. The new
[independent input audit](diverse_inputs_checked.json) verifies twenty separated
heavy-atom geometries at 0.5 Angstrom for each of twelve graph/configuration
sets, preserving full stereochemistry and using nonreflective graph-matched
alignments. These 240 starts are inputs, not 240 optimized minima or free energies.

Prepared structures, frame identities, inherited DFT reference anchors and
their hashes are in
`results/polymer_models/campaigns/2026-10-25_diverse_conformer_inputs/`.
Only ended source searches are read; completed and timeout-checkpoint sources
remain explicitly distinct. Search hit counts never become Boltzmann weights.
Both audits bind their source geometries and the RDKit version. Raw searches
remain at the source locations recorded in the individual plans.

The [40-case empty-backbone geometry campaign](../../polymer_models/campaigns/2026-10-25_smooth_empty_geometry_batch/README.md)
is a separate experimental-model calculation. Further occupied-graph resource
pilots, DFT energy re-ranking, duplicate-basin treatment and native ensemble
thermochemistry remain necessary before declaring sampling convergence.
