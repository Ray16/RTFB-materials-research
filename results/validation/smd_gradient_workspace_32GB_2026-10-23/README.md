# Bounded-gradient recovery: exact 115-atom point

The [safe-claim point campaign](pilot_occ_2_safe_claim/plan.json) recovers only
the geometry-level gradient and sampled derivative checks after the preserved
single-point calculations. [The final receipt](pilot_occ_2_safe_claim/gradient_qualification.json)
records four energy-difference checks: two Cartesian components at 0.0005 and
0.001 bohr. The maximum error is 9.70e-8 hartree/bohr, below the unchanged 1e-6
gate. Peak GPU memory was 20,304 MiB on a 32 GB V100. The gradient maximum differs
from the earlier H100 point by 1.68e-8 hartree/bohr.

The [independent review](independent_review.json) reconstructs all four
derivatives from eight displaced SCF energies, verifies each displaced geometry,
and checks all backing hashes. Its
[CPU checker](../../../scripts/validation/check_polymer_gradient_recovery.py)
does not call the quantum producer.

The MeCN and gas-phase hybrid-SCF HDF checkpoints were independently checked
against saved molecule/state identities and the previous point energies. Their
largest energy difference is below 8e-9 eV. They were not recomputed or promoted
to an optimized state. The recovered gradient likewise supplies no minimum,
native Hessian, free energy or candidate admission.

[Small spin controls](controls/1/status.json) and the
[doublet control](controls/2/status.json) establish implementation equivalence;
the separate large-point check establishes suitability only at this geometry.
The underlying energy solver remains direct LU. The adapter uses the installed
bounded gradient contractions without changing its energy-storage flag.

The [source-pinned producer](../../../scripts/validation/qualify_polymer_gradient_workspace.py)
records inputs, source hashes, actual ownership, device/CPU preflight and raw
outputs. Raw derivative points and SCF evidence remain under
`calcs/validation/smd_gradient_workspace_32GB_2026-10-23`.
[preparation_supersession.json](preparation_supersession.json) documents an
unlaunched preparation replaced before submission to strengthen duplicate-run
protection. No failed calculation was erased or silently relabeled.

The next [full geometry/energy branch](../../polymer_models/campaigns/2026-10-23_crest_refinement/baseline__design_pmdi_ethyl__occ_2_workspace_32GB/README.md)
has started separately. Native thermochemistry, broader sampling convergence
and the other occupied components remain pending.
