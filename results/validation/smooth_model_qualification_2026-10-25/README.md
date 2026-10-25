# Experimental smooth solvent-area model controls

[Combined qualification](qualification.json) binds a distinct experimental model:
PySCF switching/Gaussian solvent-accessible-area CDS, 2,030 points per sphere,
with the unchanged 590-point SMD electrostatic surface in MeCN. The portable
method identity is `2bfdff2bfdd7d8728d9398ddaf9eeed7f8d3171c7cbceb71ed1083070945a280`.
It does not replace or relabel the primary model's numerical records.

Checked evidence comprises [component scans and quadrature](../smooth_CDS_2026-10-25/pilot/summary.json),
[all-component gradients](../smooth_CDS_2026-10-25/full_gradient_panel/summary.json),
and four complete neutral quantum controls, independently reviewed in the
[C/H/O](../smooth_quantum_2026-10-25/independent_review.json) and
[N/Cl](../smooth_quantum_NCl_2026-10-25/independent_review.json) panels.
Full quantum finite differences test twelve force components and two Hessian
columns at both steps, with independently reconstructed modes/nuclear sums.
The controls span formaldehyde, methyl radical, ammonia and chloromethyl radical.

The [CDS curvature figure](../../figures/validation/CDS_surface_smoothness.pdf)
and its [data manifest](../../figures/validation/CDS_surface_smoothness.json)
compare the legacy contact cusp with the separate smooth component. Raw arrays,
SCF checkpoints and force-displacement evaluations remain under
`calcs/validation/smooth_quantum_2026-10-25/` and
`calcs/validation/smooth_quantum_NCl_2026-10-25/`.
Their hashes and source snapshots are recorded by the independent reviews.

A [91-atom geometry pilot](../../reviews/2026-10-25/smooth_geometry_pilot/summary.json)
also passes graph/stereochemistry and the unchanged final-force checks. That
resource/geometry result supplies no active single-point energy, native Hessian
or G. Charged full-quantum coverage, polymer free energies, predictive solvent
accuracy and conformer convergence remain unqualified. Every alternate-model
reaction/occupancy cycle must use consistently recalculated components.
