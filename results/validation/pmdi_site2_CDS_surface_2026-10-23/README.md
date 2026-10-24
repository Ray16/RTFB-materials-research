# Site-2 PMDI force discrepancy: isolating the CDS term

The 115-atom optimizer returned a geometry and an active single point, but its
final total-force check failed. Both sampled 0.0005-bohr derivatives agree with
the analytic force within 1.11e-8 Ha/bohr; at 0.001 bohr, the discrepancies are
about 1.92e-5 Ha/bohr, above the unchanged 1e-6 criterion. That workflow remains
unaccepted and cannot supply the candidate's polymer cycle or native G.

The [checked component diagnostic](checked/summary.json) isolates the legacy
SMD cavity/dispersion/solvent-structure (CDS) term at the exact saved geometry,
using seven finite steps. It executes the installed GPU4PySCF package's CPU-only
Fortran CDS library through its original ctypes wrapper, without importing
GPU4PySCF/CuPy or running an SCF. The PySCF package's separate CDS library gives
identical energies at all 29 sampled coordinates. The wrapper/library source
paths and hashes, input geometry and stored quantum points are retained.

The CDS coarse-minus-fine derivative changes account for the sampled total-force
discrepancies. After isolating this term, the largest remaining difference is
1.11e-8 Ha/bohr. This is an attribution of the stored discrepancy; the remainder
is not a newly computed full quantum gradient. The two finest CDS derivatives
agree within 1.05e-12 Ha/bohr. No tolerance or production state was changed.

[Independent reconstruction](independent_review.json) uses only standard Python
and checks all 29 coordinate/energy samples, 14 CDS derivative pairs and the four
original full-energy comparisons. Reproduce it with:

```bash
python3 -I -S scripts/validation/review_cds_surface.py \
  results/validation/pmdi_site2_CDS_surface_2026-10-23/checked \
  --out NEW_REVIEW_FILE.json
```

The [two-panel figure](../../figures/validation/pmdi_site2_CDS_surface.pdf) shows
the step dependence and component isolation. Its [manifest](../../figures/validation/pmdi_site2_CDS_surface.json)
binds the input/review and displayed exports. The initial layout is preserved
under `calcs/figures/pmdi_site2_CDS_surface_2026-10-23/first_layout/`; version-2
exports supply the stable aliases. Color layout has been visually inspected.

The [separate fine-step quantum diagnostic](../pmdi_site2_fine_force_2026-10-23/README.md)
tests the full force with three smaller displacements and repeated baselines.
Even a passed diagnostic will require explicit protocol review and subsequent
native physical checks; it cannot erase the original failure or qualify G.

Producers: [component calculation](../../../scripts/validation/diagnose_cds_surface.py),
[independent replay](../../../scripts/validation/review_cds_surface.py) and
[plotter](../../../scripts/plotting/validation/plot_cds_surface_diagnostic.py).
The first parser attempt stopped before any energy call because another solvent
table entry was an expression; its source/failure note is preserved separately.
No quantum job or source state was affected.

The subsequent full quantum fine-step diagnostic passed all declared criteria. Its separately named electronic qualification preserves the failed coarse result and changes only the sampled force-validation procedure; it does not claim a new geometry, hybrid-energy calculation, native G or synthesis eligibility.
