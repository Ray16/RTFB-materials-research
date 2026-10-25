# Legacy CDS contact diagnostic

[Measured component scan](alternating_conf005.json) · [Reproduction script](../../../scripts/validation/diagnose_CDS_contact.py)

A failed alternating-backbone geometry has the implicated carbon pair at 4.199893 Å. Moving one carbon radially across 4.2 Å produces a jump of 3.4365e-4 hartree/bohr between one-sided slopes of the legacy CDS energy. The centered curvature grows approximately inversely with displacement rather than approaching a finite value. This supports a non-differentiable contact in this component, so repeated optimization or a smaller Hessian step alone cannot establish harmonic thermochemistry at that contact.

The scan evaluates only the stored legacy CDS component on a CPU. It neither repeats full DFT nor qualifies a replacement solvent treatment, a geometry, a Hessian or G. The failed structures remain excluded. An independently validated treatment of these contact minima is needed before claiming broader conformer convergence; accepted sampled minima are not silently replaced or relabeled.

A separate [final graph and stereochemistry audit](geometry_review.json) finds no structural issues. This does not establish stationarity or repair the nonsmooth energy component.
