# Saved site-2 electronic state: separately qualified force refinement

The original workspace workflow remains failed and unchanged. This new capsule
qualifies its saved geometry and active electronic energy under the explicit
`CDS_explained_step_converged_total_force_v1` validation policy. It does not claim
that the original 0.001-bohr check passed, and performs no new optimization or
hybrid-energy calculation.

The policy requires quantitative CDS attribution of the coarse discrepancy,
matching exact geometry/source hashes, and full quantum checks at 0.0005,
0.0001, 0.00005 and 0.00001 bohr for both original components. The latter three
steps are independently replayed, with three baseline repeats. The absolute
force-error bound remains 1e-6 Ha/bohr; repeatability and finest-step agreement
are additional requirements. An unexplained discrepancy, missing component,
failed refinement or changed geometry is refused.

[qualification_review.json](qualification_review.json) retains the original
failed check and both independent diagnostic reviews. The largest new full-force
error is 2.42e-7 Ha/bohr, the largest finest-pair difference 2.14e-7 Ha/bohr,
and baseline energy spread 1.82e-12 Ha. The checked stationary force maximum
is approximately 5.57e-5 Ha/bohr. These are sampled derivative checks, not a
full-Hessian or minimum certificate.

The [new electronic receipt](energy/done/baseline__design_pmdi_ethyl__occ_2__pose0.json)
retains the exact original active energy and its backing bytes. [plan.json](plan.json)
records read-only source storage and the changed validation procedure. No source
record is rewritten. Subsequent native work uses separate output storage and
retains all original Hessian, finite-gradient, imaginary-mode and partition gates.
No G or ensemble convergence is created by this capsule.

Producer: [qualify_saved_polymer_energy.py](../../../../../scripts/pipeline/qualify_saved_polymer_energy.py).
The [component audit](../../../../validation/pmdi_site2_CDS_surface_2026-10-23/README.md)
and [fine-step quantum diagnostic](../../../../validation/pmdi_site2_fine_force_2026-10-23/README.md)
retain raw paths, source hashes and independent replay commands. The ordinary
candidate observer handles the next native-thermal stage; full-axis eligibility
and synthesis admission remain false.
