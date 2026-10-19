# First independently checked catholyte parent addition

Exact free N-ethylphenothiazine (EPT), first one-electron oxidation / +1-to-0
reduction couple: computed **+0.537375 V vs Fc**, experimental **+0.278 V vs Fc**,
residual **+0.259375 V** before display rounding. The comparison is positive-side
parent validation, not a measured property of any grafted structure.

The source is Attanayake et al., Chem. Mater. 2019, DOI
10.1021/acs.chemmater.8b04770, Table 1, p.4356. CV conditions are 1 mM analyte,
0.5 M TEATFSI/acetonitrile, 26 C, 10 mV/s; the authors calibrate using ferrocene.
The model uses 298.15 K thermochemistry and continuum MeCN; supporting salt is not
modeled explicitly. The second wave is excluded from strict scoring because the
source does not establish chemical reversibility.

Both endpoint native records and the live matched Fc reference are independently
read. Both relevant alternate-spin optimizations converge; the neutral triplet
and cation quartet are 3.09384 and 3.60706 eV above the selected surfaces. These
are electronic comparisons with stated thermal uncertainty, not a full spin or
conformer ensemble. Six independent arithmetic/admission comparisons pass with
stable inputs, including explicit experimental CSV hashes.

The single point overestimates experiment. It cannot determine a family error
bar or justify an offset fit. The earlier seven-point native benchmark is a
historical selected-spin snapshot: two PMDI second waves now await the existing
DFT spin gate. Their targeted calculations are under
`results/validation/benchmark_imide_spin_2026-10-19/`. Aggregate accuracy is not
updated by silently retaining those unresolved points.

The new reusable review is `scripts/validation/benchmark/review_matched_parents.py`;
it writes a fresh scoped directory and performs no quantum calculations or
production-table promotion. Pending/source-excluded comparisons stay blank.
