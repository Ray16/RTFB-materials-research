# Experimental DIPEA salt pKa and base-transfer sensitivity

[Source-bound assessment](summary.json) · [Focused arithmetic and source tests](tests.xml) ·
[Prior 18.6 source-origin audit](../dipea_reference_origin_2026-10-23_checked/README.md)

Cook et al. report an NMR titration/competition **estimate of 18.7** in MeCN
for diisopropylethylammonium **hexafluorophosphate**:
[Chemical Science, 2024, DOI 10.1039/D3SC02611F](https://doi.org/10.1039/D3SC02611F).
The exact primary article and measurement statement were checked from the
Europe PMC full-text XML, whose byte hash is in the assessment. The raw XML is
retained under `calcs/literature/2026-10-26_dipea_experiment/` and is CC BY-NC
3.0, distinct from the repository's code license. We do not infer the experiment's
exact temperature, ionic activities or salt dissociation from the main article.

For a **conditional cross-source comparison**, combine the reported 18.7 with
the existing DBUH pKa 24.31 and HCl pKa 10.30 at 298.15 K. The relative base
transfer `DIPEA + DBUH+ -> DIPEAH+ + DBU` is +31.80 kJ/mol in the stored DFT
components and +32.02 kJ/mol from these pKa values. This comparison removes the
common HCl/chloride terms, but the experimental salt and reference scale are
not demonstrated to have identical activity and ion-pairing conventions; the
0.22 kJ/mol difference is not an uncertainty estimate.

Under the same cross-source assumptions, the pKa-implied separated-ion
`DIPEA + HCl -> DIPEAH+ + Cl-` capture term is −47.95 kJ/mol, versus −70.98
kJ/mol from DFT. Their difference is −23.03 kJ/mol in the DFT-minus-reference
direction. HCl cancels in the complete all-DFT grafting cycle. Replacing only
the capture term changes shared-error cancellation and produces **hybrid
sensitivities**, not validated grafting predictions. For the two BQ cases those
sensitivities are approximately +9.69/+8.08 kJ/mol, while the original all-DFT
DIPEA sums remain −13.34/−14.94 kJ/mol. No reaction result or route decision
is silently replaced, and the observed agreement of relative base transfer
does not uniquely identify the cause of the capture discrepancy.

The older 18.6 value remains labeled as a PCM estimate in its author-posted
source. The newer experimental salt estimate is separate evidence; it does not
retroactively make that older number experimental. The NMR estimate supplies
an informative benchmark, while chloride speciation, homoconjugation,
temperature/activity matching and full grafting accuracy remain unresolved.
