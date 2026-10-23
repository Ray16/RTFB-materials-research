# Fixed-geometry nuclear rotor-cutoff sensitivity

[summary.json](summary.json) propagates the archived projected nuclear partitions
at 25, 50 and 100 cm^-1 through all 39 benchmark endpoints and 87 completed
grafting reaction models (67 HCl-reference and 20 base-assisted). Geometry,
electronic energy, Hessian, spin factor and standard-state convention are fixed.
The primary 50 cm^-1 results are preserved. No new SCF calculation is performed.

The largest endpoint span is 25.43 mV; the largest displacement from the primary
setting is 23.01 mV. Aggregate MAE is 0.2449/0.2470/0.2473 V at the three cutoffs.
The Fc reference itself spans 2.82 mV. The largest reaction span is 4.46 kJ/mol;
none of the 87 nominal reaction signs changes. A specified base/byproduct term
remains constant across its five candidates at each cutoff. New variant reaction
sums use the exact SI conversion; original source arithmetic is not overwritten.

The [independent full-state reconstruction](independent_review.json) checks 117
potential evaluations and 261 reaction sums, rather than reusing the producer's
relative-shift formula. Maximum differences are 1.15e-11 V and 8.66e-10 kJ/mol.
Both scripts verify every source-package byte against its inventory. The local
raw-state source is `calcs/releases/2026-10-23_scientific_review_v5/`, whose archive
and exact scope are indexed in [the package guide](../../releases/primary_review/README.md).

Producer: [thermal_cutoff_sensitivity.py](../../../scripts/validation/thermal_cutoff_sensitivity.py).
Independent checker: [review_thermal_cutoff_sensitivity.py](../../../scripts/validation/review_thermal_cutoff_sensitivity.py).
Plotter: [plot_thermal_cutoff_sensitivity.py](../../../scripts/plotting/validation/plot_thermal_cutoff_sensitivity.py).
The stable [PDF figure](../../figures/validation/thermal_cutoff_sensitivity.pdf)
and [manifest](../../figures/validation/thermal_cutoff_sensitivity.json) retain
source hashes and displayed data. Version-2 immutable exports are also retained;
the original layout is preserved under `calcs/figures/thermal_cutoff_sensitivity_2026-10-23/first_layout/`.
The figure was inspected in color and grayscale and at SI size.

These ranges quantify a settings sensitivity. They do not establish candidate
prediction intervals, conformer convergence, ion speciation, selective attachment
or experimental accuracy. No all-axis candidate is admitted by this check.
