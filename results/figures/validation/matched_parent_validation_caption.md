# Matched-parent redox-potential benchmark

Panel a compares the measured and calculated potentials (V versus Fc/Fc+ in
MeCN). The diagonal is agreement, with no fitted shift. Panel b shows signed
calculated-minus-measured errors for each parent; symbols distinguish the first
and second redox steps. Twelve tier-A potentials from eight exact free parents
and four studies give MAE 0.20904 V and RMSE 0.24572 V. DMeOEPT's second step is
0.655 V measured and 1.202 V calculated, an error of +0.547 V.

Statistics reproduce the snapshot's three-decimal calculated-potential
serialization. Wave spacings are derived quantities and are excluded from the
12-event error statistics. Unadmitted second waves remain absent. States retain
the spin-confidence qualifications in the source comparison table. Measurements
from the same parent are correlated; this small set does not supply a universal
uncertainty bound. Supporting electrolytes are not modeled explicitly and PMDI
salt assignment remains ambiguous. The result does not establish accuracy for
exact grafts, polymers, reorganization energies or electron-transfer rates.

Data and experimental source DOIs:
`results/reviews/2026-10-19/full_matched_parent_final/comparisons.csv`.
The filtered actual values are in `matched_parent_validation.csv`; PDF/SVG are
vector exports and PNG is 600 dpi. Regenerate with
`python scripts/plotting/validation/plot_matched_parent_validation.py`.
