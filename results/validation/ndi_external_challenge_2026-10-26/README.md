# NDI external experimental challenge: source and scope

[Source-bound record](summary.json) · [Independent checks](independent_check.json) ·
[Tests](tests_final.xml) · [Frozen 39-endpoint benchmark](../publication_benchmark/README.md)

[Nowacki and coauthors (2025)](https://doi.org/10.3390/ijms26104807) report
cyclic voltammetry of an NDI–styrene molecule in acetonitrile with 0.1 M
tetrabutylammonium perchlorate, calibrated against external Fc/Fc⁺. The
reported cathodic/anodic peaks are −1.040/−0.880 V for the first step and
−1.440/−1.350 V for the second. Their arithmetic midpoints are **−0.960 and
−1.395 V**, with peak separations of 0.160 and 0.090 V. These midpoints are
derived proxies: the article describes quasi-reversible processes and does not
report them as formal potentials. The figure caption unambiguously specifies
100 mV/s for the polymer trace, not the monomer; this record does not assign
that scan rate to the monomer.

The proposed exact free-molecule graph has formula C30H18N2O4 and calculated
exact mass 470.1267 Da, consistent with the paper. Its two para-vinylphenyl
imide attachments match the authors' synthesis and Figure 1, which was inspected
independently; Figure 5A was also inspected for the solution voltammogram. The
graph is absent from the checked project molecular manifest and the 837-design
search space. That is a project-local identity check, not a novelty claim.

[Ten checked starting geometries](../ndi_external_control_seeds_2026-10-26/README.md)
are now prepared; they are force-field inputs and do not change the experimental
or benchmark status below.

This is a **free-molecule method control**, not a proposed pre-polymerization
synthesis route or a model of a Merrifield graft. Its charge/spin states are
proposals, and no computed potentials are available. It is not added to the
frozen 39-endpoint panel or used in its accuracy metrics. To assess it further,
the exact molecule requires checked geometries, alternative-spin assessment,
native thermochemistry and a method comparison; the peak midpoint's suitability
as an equilibrium target must also be considered. It cannot validate a different
NDI derivative or polymer environment by structural analogy alone.

The XML and independently inspected source images are retained under
`calcs/literature/2026-10-26_ndi_source/`; the primary article is CC BY 4.0.
The record stores their hashes. A superficially relevant [Erten et al. (2005)](https://doi.org/10.1016/j.jphotochem.2005.04.038)
NDI electrochemistry study used **chloroform** for its cyclic voltammetry;
its acetonitrile fluorescence measurements do not make those redox potentials
MeCN endpoints.
