# Independent DIPEA protonation measurements in acetonitrile

[Computed comparison](summary.json) · [Independent source review](independent_review.json) ·
[Focused tests](tests.xml) · [CPU regression and separately resolved skips](full_regression_summary.json) ·
[Earlier NMR assessment](../dipea_experimental_anchor_2026-10-26/README.md)

Barišić et al. report **log Kᴴ = 18.64(1)** for DIPEA in MeCN at
25.0 ± 0.1 °C from a spectrophotometric titration using bromocresol green.
The value is printed below Figure S51 on page 43 of their supporting
information; the main paper describes the indicator constants used in the fit.
[RSC Advances, 2021, DOI 10.1039/D1RA04738H](https://doi.org/10.1039/D1RA04738H).
Figure S41 concerns DMSO and must not be used as the MeCN value. The raw
article XML, supporting-information PDF and independently rendered page remain
on shared storage in `calcs/literature/2026-10-26_dipea_experiment/`.

The separate 2024 NMR competition estimate for DIPEAH hexafluorophosphate is
18.7 in MeCN. Its supporting information specifies MeCN-d3, an initial 0.02 M
salt solution and a triethylamine reference of 18.81. Its Table S3 contains
negative entries under a Ka-ratio heading and reports the logarithm of the
absolute ratio. We retain the authors' stated estimate and do not infer a new
signed equilibrium fit from that table. These measurements are not combined
into a pooled error bar.

At 298.15 K, using 18.64 with the separately tabulated HCl and DBU references
gives −47.60 kJ/mol for DIPEA/HCl capture and +32.36 kJ/mol for the HCl-free
DIPEA/DBUH⁺ proton-transfer comparison. The corresponding stored DFT values
are −70.98 and +31.80 kJ/mol. These are conditional cross-source comparisons:
indicator calibration, ionic activities, chloride association and the acid/base
scales are not established to be identical. The small relative-base discrepancy
is not a prediction interval or a unique diagnosis of any individual species.

HCl cancels in the full all-DFT grafting cycle. Substituting only the inferred
capture term changes error cancellation; the BQ values of approximately
+10.03/+8.43 kJ/mol are hybrid sensitivities, not validated grafting predictions.
The all-DFT values and candidate decisions remain unchanged. Missing precursor
values remain missing. Selectivity, yield, chloride speciation and the full
reaction free energy remain unvalidated.
