# Frozen local motif interactions

[Current sums](summary.json) · [Construction audit](../local_motif_construction/summary.json)
· [Parent conformers](../../../figures/polymer_models/spacing_conformers.pdf)
· [Prepared inputs](../../../polymer_models/campaigns/2026-10-26_local_motif_fragments/manifest.json)
· [Exact cut/cap construction figure](../../../figures/polymer_models/local_motif_construction.pdf)

The first complete checked triplet, adjacent-site parent pose 0 with 1.09 Å caps,
gives **−37.16 kJ/mol in gas and −24.76 kJ/mol in MeCN** for
`E(AB) − E(A) − E(B)`. Each state has its own converged single point and
independent SCF-checkpoint, geometry, charge/spin, runtime and resource review.
All three use ωB97M-V/def2-TZVPD/VV10 at the same frozen coordinates and medium.
No standard-state term is added to these electronic energies.

This model retains each phthalimide motif and its primary benzyl carbon; replacing
the bond to the polymer aryl group with H makes an N-methylphthalimide fragment.
Original motif atoms and hydrogens retain their exact parent positions. The two
fragments preserve their relative pose, with no covalent bond between them.
Full-precision XYZ files define calculation coordinates; MOL files provide bond
topology and rounded drawing coordinates. The independent audit checks both.

The two chemical graphs, one capped motif and a disconnected pair, occur in
36 geometry-specific inputs: four accepted close/separated parent poses,
three predeclared cap lengths (1.07/1.09/1.11 Å), and states A/B/AB per model.
These form twelve interaction triplets. All 24 monomer inputs now have checked
electronic records, including the original two references. The remaining
11 pair inputs have a separate dependency-gated campaign. A far-pose pair
resource sentinel precedes expansion of that pair class. Missing values stay
blank; the construction inputs are method controls and do not expand the
enumerated candidate count.

The gas and MeCN sums quantify the frozen methyl-capped local electronic
interaction, including caps and continuum response. Basis-superposition
correction and broader basis checks remain pending. These values do not supply
thermal association free energies, pure π binding, a decomposition of the
whole-chain grafting cycle, electron-transfer coupling or conformer populations.
They must not be subtracted from polymer ΔΔG to assign a backbone free energy.

Raw calculations are under `calcs/local_motif_interactions/`, in the exact
namespaces in each plan. Pilot references are in
`results/polymer_models/campaigns/2026-10-26_local_motif_energy_pilots/` and the
separate first-pair B campaign. All use the shared GPU gate and disjoint CPU
sets. The measured 19/38-atom frozen-energy class is the scope of the T4 tests;
large-polymer geometry and Hessian suitability remain separately controlled.
