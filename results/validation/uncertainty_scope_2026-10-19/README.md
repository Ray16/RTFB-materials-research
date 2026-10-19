# Computational uncertainty scope and sparse-evidence guards

The lambda scale is a descriptive spread between two computational models,
not an experimentally calibrated solution-rate error. The current comparison
contains 442 rows; 330 lie in the declared 0–2 eV range. Their MAD-derived spread
is 0.209178 eV and median residual is -0.100120 eV. The source CSV does not carry
per-record convergence, binding or conformer-QC columns. None of those qualities
is inferred from a value-range filter.

The scorecard now exposes the computational-only status, residual count,
minimum comparison count, spread, median residual and configured prior floor.
A single residual cannot estimate a spread. Zero/small computational spread
cannot tighten the configured 0.10-eV prior. Adequate spread may widen it; the
current 0.209-eV value remains unchanged. The median is diagnostic and is never
used as a fitted correction. Three regression cases cover a single observation,
identical discrepancies and a larger spread. Both primary and independent table
calculators apply the same stated minimum/floor through separate code.

The [primary CV condition check](primary_condition_review.json) rechecks the
[2019 source](https://xray.uky.edu/people/parkin/papers/426.pdf), methods p.4355:
1 mM active material, 0.5 M TEATFSI in MeCN, 26 °C, and separate 5 mM ferrocene
calibration before each measurement. The pure-solvent bare-ion model does not
explicitly represent that electrolyte. This is a condition mismatch, not proof
that ion pairing causes the model error. Global reference offsets cannot remove
the charging-gap discrepancy. No experimental label or canonical prediction is
changed by this source review.

The full suite passes 544 tests. New property tables and candidate exports are
reconstructed after the reporting-contract change; publication requires their
independent audit and unchanged input inventories. Physical accuracy and solution
kinetics remain separate, unvalidated scopes.
