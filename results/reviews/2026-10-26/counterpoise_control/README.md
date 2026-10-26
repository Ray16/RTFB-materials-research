# Gas counterpoise control and prepared motif case

[Independent CPU/GPU review](independent_review.json) · [Focused tests](tests_v2.xml) ·
[Water-dimer control](../../../validation/local_counterpoise_2026-10-26/water_dimer/plan.json) ·
[Prepared close-motif case](../../../validation/local_counterpoise_2026-10-26_motif_v2/occ_2_3__pose0__cap109/plan.json) ·
[Independent geometry review](prepared_motif_geometry_v2.json)

Five gas-phase water-dimer terms were computed on the same full-complex
quadrature and auxiliary basis. All CPU terms used the 226-function full
auxiliary space. The independent review checks the GPU and CPU checkpoints,
exact nuclei/electrons/coordinates, common grid, orbital occupancy, variational
checks, runtime ownership and source hashes. The maximum CPU/GPU total-energy
difference is 3.47×10⁻⁸ hartree; the counterpoise interaction differs by
8.14×10⁻⁵ kJ/mol. This tests numerical implementation on the stated control;
it does not establish experimental interaction accuracy or a solvent correction.

The earlier CPU attempts are preserved as unqualified. In those attempts,
PySCF's density-fitting build replaced the intended full auxiliary molecule
for isolated fragments. A scoped process-local factory now supplies the full
complex auxiliary space during actual `DF.build()`, including its range-separated
branch, and restores the installed factory afterwards. Focused tests exercise
the real build path and source restoration; installed packages are unchanged.

The 38-atom close phthalimide control has been prepared from the exact reviewed
fragment coordinates. Its independent geometry check verifies 168 electrons in
the complex, 84 for either ghost monomer and 1,026 full/ghost AO functions.
No scientific job has started for this case. Its guarded launcher requires a
32-GB-class GPU, checked both before launch and inside the reserved worker, and
retains the ordinary reservation, CPU, host-memory and residency checks. The
earlier unsubmitted plan is preserved with its source and a supersession record.
The exact plan and the independently checked small control are
required before launch; unavailable devices do not consume a scientific attempt.

Raw control calculations remain in
`calcs/validation/local_counterpoise_2026-10-26/water_dimer/` and
`calcs/validation/2026-10-26_counterpoise/`. The prepared point will write to
`calcs/local_motif_counterpoise/`. No counterpoise adjustment is applied to
MeCN, thermal free energy, the polymer cycle or a candidate-selection decision.
