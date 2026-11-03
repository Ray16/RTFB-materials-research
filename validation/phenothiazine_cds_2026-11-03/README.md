# Phenothiazine CDS derivative accuracy and native validation

The two neutral full Hessians failed the unchanged 1e-6 Ha/Bohr-squared
symmetry limit. CPU reconstruction identifies the legacy SMD CDS numerical
Hessian as the source of the antisymmetry. The saved full matrices have maxima
6.543e-6 and 6.562e-6; removing the independently recomputed default CDS term
leaves residual asymmetries of 5.13e-13 and 5.56e-13. The original failed
matrices and their calculation records remain intact.

The [minus](minus_component.json) and [plus](plus_component.json) component
studies use the installed legacy function and its Fortran library, without
importing CuPy or starting CUDA. Both exact native centers are retained.
Steps of 1e-5, 1e-6 and 1e-7 Bohr give component/reassembled asymmetries below
1e-9. Energy differences check all 117 CDS gradient components at two steps;
the maximum discrepancy is below 4.61e-10 Ha/Bohr. These component studies
qualify a [targeted numerical step control](qualification_02.json), not a complete
Hessian, minimum, G or redox potential.

Three [fresh full native controls](prepared_manifest.json) use a 1e-6-Bohr CDS
gradient-difference step with the same energy and gradient functions. They
recompute the full DFT Hessian; the diagnostic reassembled matrices are not
substituted into production results. The worker records the actual CDS
evaluation, center, helper hash and numerical policy. The original 32-GB-class
resource envelope and all stationarity, symmetry, minimum, spin, two-step
derivative and sulfur-column gates remain mandatory. [Archive preflight](archive_preflight.json)
passed for all three private namespaces. The [current record](current.json)
distinguishes queued work from scientific starts and qualified properties.

The [original oxidized native state](original_oxidized_native_review.json)
independently passes all 21 actual SCF/frame/spin checks, standard and sulfur
derivatives, masses, projection and nuclear terms. Its minimum frequency is
26.0728 cm-1 and maximum sampled-column relative error is 2.3382e-5. This is
one qualified component at the default 1e-4-Bohr CDS step; it does not supply
the missing neutral G or an oxidation potential.

The independent partition reader initially looked up the current refactored
qRRHO wrapper instead of the producer's declared archive. The archive retains
the exact source hash required by its certificate. The reviewer now binds that
explicit archive, preserving all certificate, library, artifact and numerical
checks. It reads the independently checked spin term from the full projected
record. No quantum output was altered to resolve these reader issues.

Use `python -m redox.validation.cds_hessian_component --help` for the CPU
diagnostic and `python -m redox.workflows.route_cds_native --help` to reconstruct
the component qualification. The maintained native workflow accepts the explicit
`--cds-component-qualification` argument for this scoped control. Numerical
policies remain attached to reviewed G values; any subsequent reference or
cycle comparison must retain them. The [130 focused tests](final_checks_02.xml) and
the actual oxidized-state reconstruction are separate validation evidence.

The queued namespaces preserve the original qualification receipt and its
[as-used qualifier source](qualification_producer_v1.py). The second
qualification adds explicit nonfinite-array and gradient rejection and
reconstructs identical scientific results. This review does not alter the
queued numerical protocol or overwrite its frozen sources.
