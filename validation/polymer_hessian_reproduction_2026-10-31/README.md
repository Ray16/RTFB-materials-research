# Full corrected polymer-Hessian reproduction

Three complete corrected Cartesian Hessians for the same retained 119-atom
phthalimide control pass independent physical reconstruction and every
prospectively fixed pairwise numerical comparison. The matrices originate from
a fresh V100 calculation, a V100 calculation reusing completed solute components
while recomputing solvent H1, and the matched H100 component-reuse calculation.
Fresh displaced-gradient and nuclear/mode validation was completed on V100s in
separate saved-matrix namespaces. See [controls and scope](summary.json).

All sampled derivative checks pass the unchanged 1% bound; the largest relative
error is 9.45e-6. Every reconstructed internal mode is positive, with a lowest
frequency of approximately 9.59 cm-1. The 351 internal modes exclude six rigid
motions geometrically. Independent nuclear equations reconstruct all three
rotor-cutoff terms. The maximum nuclear-term difference from the retained
uncorrected reference is about 1.52e-8 eV for this structure.

The [pairwise comparison](comparison.json) passes every
[predeclared tolerance](declared_tolerances.json). Across the three calculations,
maximum differences are approximately 1.33e-9 Ha/Bohr2 for Hessian elements,
6.41e-6 cm-1 for internal frequencies and 1.67e-8 eV for nuclear terms. The
retained state, numerical protocol, actual source hashes, source response rows,
fresh SCF checkpoints, gradients and complete derivative coverage are checked.

The original matrix-producing runs stopped on wall/memory limits after retaining
the complete matrices. Their failed resource records remain intact. A completed
validation stage loads the exact source matrix, performs a fresh center solve,
checks state/geometry agreement and runs the ordinary derivative, projection,
minimum and nuclear tests. Source response records cover all 357 Cartesian rows;
no new response solve is falsely claimed by the validation stage. The source
matrix producer's component reuse is recorded separately from the validation
stage, which reuses the already complete Cartesian matrix.

This establishes numerical/physical reproduction for one large retained model.
It does not qualify a complete uninterrupted resource profile, generalize the
implementation sensitivity to every state, complete reaction cycles, assemble
new total chemical free energies or close polymer conformer convergence. The
cause of the earlier gradient discrepancy is not uniquely assigned by this
comparison. Electronic benchmark and small-molecule controls remain separate
supporting evidence.

Frozen independent reports, nuclear data and derivative receipts are sealed by
[manifest](manifest.json), with as-used producer snapshots. Large source arrays,
SCF checkpoints and raw gradients remain at the bound namespace locations.
The existing comparison script reconstructs the published differences from those
independently reviewed controls; no figure or launch receipt substitutes for the
physical review.
