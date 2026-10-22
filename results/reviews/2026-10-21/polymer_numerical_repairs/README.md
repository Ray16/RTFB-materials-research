# Polymer numerical execution checkpoint

This checkpoint verifies numerical controls, not experimental predictive accuracy.
The focused regression set passes 79 tests. Bounded solvent-Hessian contractions
reproduce frozen full singlet/radical Hessians and thermal corrections; the largest
thermal difference is 1.63e-8 eV. A same-geometry single-point replay agrees within
8.74e-11 eV before the failed alternative conformer is admitted.

[Evidence and source hashes](summary.json) retain immutable snapshots of controls
and plans. Large Hessians and calculation logs remain on shared storage at the
recorded raw locations. Optimizer, SCF, connectivity and native thermal admission
criteria remain unchanged. No absent thermal correction is replaced by zero.

Large-model native free energies and two original PMDI trimers remain pending.
The completion controller independently reconstructs matched cycles and requests
a new thermal calculation if the recovered conformer changes the electronic
minimum; it does not transfer a correction between geometries.

The current experimental benchmark is separate: see the
[validation guide](../../../validation/README.md). Repository-wide cleanup remains
deferred until the scientific workflow is reviewed.
