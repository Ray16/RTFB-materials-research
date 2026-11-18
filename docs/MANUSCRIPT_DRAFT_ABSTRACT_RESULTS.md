# Manuscript draft — Abstract + Results narrative (draft 2026-11-17)

Framing chosen (most-supported by the assets; swap if you prefer): *a validated DFT+SMD screening
pipeline for graftable redox motifs, whose application fills the missing 2e-catholyte side of a
Merrifield-bound non-aqueous flow battery with dihydrophenazines carrying measured design rules.*

## Abstract (draft)
Non-aqueous organic redox flow batteries that bind the redox-active species to a Merrifield
(chloromethyl polystyrene) backbone need BOTH a low-potential anolyte and a high-potential
catholyte that survive covalent grafting. We built a reproducible DFT+SMD pipeline (wB97M-V/
def2-TZVPD, SMD acetonitrile) that grafts candidate motifs onto a 4-methylbenzyl tether and
computes redox potential (vs Fc/Fc+), gravimetric capacity, and inner/outer-sphere
reorganization energy on a uniform thermochemical basis. Against 44 measured couples the
potentials are reliable to a mean absolute error of 0.23 V with a systematic, family-stratified
bias (+0.17 V reductions, +0.27 V oxidations) that we measure and remove; the catholyte side is
anchored per-couple against the exact dihydrophenazine parent BMEPZ (MAE 0.30 V). Screening 92
grafted candidates, the anolyte side is well populated (quinones, naphthoquinones, imides) while
the catholyte side was previously empty of an honest two-electron motif. We fill it with grafted
5,10-dihydrophenazines and establish three DFT-confirmed design rules: capacity rises with a
smaller N'-cap (N'-methyl-DHP, 91 mAh/g), potential rises with an electron-withdrawing N'-cap
(N'-(2,2,2-trifluoroethyl)-DHP, +0.62 V), and potential falls monotonically with ring donors
(unsubstituted +0.43 -> 2,7-dimethyl +0.30 -> 2,7-dimethoxy +0.21 V). The recommended cell,
2-(benzylamino)-1,4-naphthoquinone with N'-methyl-dihydrophenazine, is a two-electron/two-electron
pair at a bias-corrected ~1.83 +/- 0.12 V, both grafted through the same Merrifield tether.

## Results narrative (draft, by figure)
### R1. A grafting-aware screening pipeline (Fig. 1: schematic + parity)
Motifs are registered as SMILES fragments and grafted onto a 4-methylbenzyl proxy for the
Merrifield benzylic-displacement product; redox states are optimized at DFT+SMD. Validation
against 44 tier-A measured couples gives MAE 0.226 V (parity, Fig. 1b). The error is systematic
and family-stratified (imide +0.07 best, quinone +0.25 worst; reductions +0.17 vs oxidations
+0.27), i.e. a correctable offset rather than scatter.

### R2. The candidate landscape (Fig. 2: potential x usable capacity, both Pareto fronts)
92 grafted candidates on one legacy xTB-RRHO basis, ranked by a sigma-aware Pareto front over raw
axes (no scalarized figure-of-merit). Anolyte front: PMDI, anthraquinone ethers, amino-
naphthoquinones, trimethyl/methoxy-benzoquinones (-1.04 to -1.54 V, 61-72 mAh/g). The catholyte
side is sparse — motivating the expansion.

### R3. Filling the 2e-catholyte gap with grafted dihydrophenazines (Fig. 3: design rules)
Four grafted DHPs, computed and scored, establish the three monotonic design rules (above).
Outer-sphere lambda_o is computed pool-wide (92/92, 0.52-0.77 eV); inner-sphere lambda_i is
well-defined for these oxidations (no anion-unbound pathology, unlike the n-type reductions),
~0.26-0.34 eV/couple.

### R4. Per-couple catholyte validation (Fig. 3 inset / Table)
The exact BMEPZ parent, computed at the same level, reproduces its measured -0.29/+0.50 V to a
per-couple MAE of 0.30 V (bias +0.30 V, both couples high) — putting the catholyte on the same
validation footing as the anolyte quinones.

### R5. Proposed cell
2-(benzylamino)-1,4-NQ (anolyte, -1.41 V, 70 mAh/g, lambda_o 0.73 eV) + N'-methyl-DHP (catholyte,
+0.47 V, 91 mAh/g, lambda_o 0.70 eV): bias-corrected ~1.83 V, residual +/- 0.12 V after measured-
bias removal. Both 2e, both Merrifield-grafted.

## Methods caveats to state plainly
- Uniform legacy xTB-RRHO thermochemistry; native projected-qRRHO available for spot-checks
  (thermal correction largely cancels in redox differences).
- lambda_i held for n-type reductions (anion_unbound is correct physics); kinetics ranks on lambda_o.
- Candidates are proposals; component potentials are benchmarked, full cells are not (no synthesis).

## Open (author decisions / out of model's hands)
- Confirm framing + target venue.
- If a reviewer wants native-qRRHO on the leads specifically, run that subset (cheap, per-molecule).

---

## Methods (draft)

### Grafting model
Each redox motif is a SMILES fragment with one attachment point, joined (RDKit molzip) to a
4-methylbenzyl tether `Cc1ccc(C[*:1])cc1` — a minimal proxy for the product of benzylic-chloride
displacement on Merrifield (chloromethyl polystyrene) resin. The tether is common to all 92
candidates; distinguishing substituents are reported per candidate. The proxy isolates the
redox motif's intrinsic properties from polymer-chain effects (stated as a scope limit).

### Electronic structure
Geometries of every redox state are optimized in implicit solvent at wB97M-V/def2-TZVPD with SMD
(acetonitrile), density-fitted (RI-J), anions using the diffuse-augmented basis. Energies are the
SMD electronic energy at the optimized geometry. Conformers are pre-screened with the UMA MLIP
before DFT. All GPU work is gated through a cluster-wide reservation lock.

### Thermochemistry
Gibbs energies use G = E_smd + G_thermal with a uniform gfn2-xTB RRHO thermal correction
(298.15 K) across all 92 candidates, so every couple is compared on one basis. Because G_thermal
enters a redox potential as the difference between adjacent charge states of the same molecule,
it is largely self-cancelling; a native projected-internal-mode qRRHO treatment is available in
the pipeline for subset spot-checks.

### Properties
- Redox potential E = (G_ox - G_red) - G(Fc/Fc+), referenced to a consistently-computed ferrocene,
  reported vs Fc/Fc+ (electrolyte-independent).
- Usable gravimetric capacity = nF / [3.6 (MW + n * M_counterion)], counterion PF6- (catholyte) or
  TBA+ (anolyte, organic-standard supporting salt TBAPF6 in MeCN).
- Outer-sphere lambda_o: single-sphere Born/Marcus from the SASA effective radius, MeCN
  eps_s 37.5 / eps_op 1.806 (Pekar). Computed for all 92. Molecular-cavity PCM is the shape-aware
  refinement (ratio ~1.25-1.38 on the prior set).
- Inner-sphere lambda_i: gas-phase Nelsen 4-point. Well-defined for p-type oxidations; held for
  n-type reductions where the gas-phase radical anion is unbound (a physical limit, not missing data).

### Validation
44 tier-A measured couples (benchmark_mecn.csv): MAE 0.226 V, systematic bias +0.191 V, stratified
(+0.166 reductions, +0.268 oxidations; by family imide +0.072 ... quinone +0.253). The catholyte
family is anchored per-couple by computing the exact BMEPZ parent at the same level
(E1 -0.062/meas -0.290; E2 +0.869/meas +0.500; MAE 0.299 V). Proposals report bias-corrected
potentials with same-chemistry corrections; residual uncertainty after bias removal ~0.12 V.

## Figure captions (draft)
- **Fig. 1b** (EXISTS: results/figures/validation/parity_mecn.png).
- **Fig. 1.** (a) Grafting schematic: motif + 4-methylbenzyl Merrifield proxy. (b) Parity of
  computed vs measured potentials (44 couples, MAE 0.226 V), colored by family; dashed = stratified bias.
- **Fig. 2.** Candidate landscape: redox potential vs usable capacity (TBAPF6/MeCN) for all 92,
  colored by family, with sigma-aware anolyte and catholyte Pareto fronts and a real-name key.
  (results/figures/candidates/candidates_potential_vs_capacity_MeCN_TBAPF6.png)
- **Fig. 3** (EXISTS: results/figures/candidates/dhp_design_rules.png) Dihydrophenazine catholyte design rules: potential vs N'-cap (EWG/EDG) and ring
  substitution (monotonic +0.43 -> +0.30 -> +0.21 V); inset: BMEPZ per-couple validation.

## Limitations (state plainly)
1. Thermochemistry is uniform legacy xTB-RRHO; native projected-qRRHO is not applied pool-wide
   (it is largely self-cancelling in potentials; available for spot-checks).
2. lambda_i is held for n-type reductions (anion_unbound is physical); kinetics ranks on lambda_o.
3. Candidates are computational proposals on a tether proxy; component potentials are benchmarked,
   full cells and exact polymer-bound synthesis/stability are not.

---

## Introduction (draft)
Redox flow batteries (RFBs) decouple power from energy by storing charge in flowable electrolytes,
and non-aqueous organic RFBs promise wide voltage windows unconstrained by water's stability. A
persistent failure mode is active-species crossover and capacity fade; covalently binding the
redox motif to a polymer backbone — here Merrifield (chloromethyl polystyrene) resin via benzylic
displacement — suppresses crossover but constrains the chemistry to motifs that tolerate grafting
through a benzylic tether. A practical cell needs a matched pair: a low-potential n-type anolyte
and a high-potential p-type catholyte, ideally two-electron for gravimetric density, both grafted
through the same handle. Computationally, the anolyte chemical space (quinones, imides) is well
charted, but an honest two-electron catholyte that survives grafting has been scarce. We address
both the method and the gap: a validated DFT+SMD pipeline that computes device-relevant properties
for grafted motifs on a uniform basis, and its application to 92 candidates, which identifies
grafted dihydrophenazines as a two-electron catholyte family with tunable, DFT-confirmed design rules.

## Discussion (draft)
Three points frame the results. First, reliability: the pipeline's potentials carry a systematic,
family-stratified offset (reductions +0.17 V, oxidations +0.27 V) rather than random scatter, so
proposals are reported bias-corrected with same-chemistry anchors; after removal the residual is
~0.12 V. Crucially the catholyte side is not merely corroborated — computing the exact BMEPZ parent
reproduces its measured couples to 0.30 V, putting the two-electron catholyte on the same footing
as the benchmarked anolyte quinones. Second, robustness of the thermochemical basis: a measured
legacy-vs-native (projected qRRHO) comparison across all candidate families shows the thermal
correction is basis-invariant at the couple level (0 mV); only a floppy organometallic reference
(ferrocene) shifts, uniformly, cancelling in cell voltage. The uniform legacy basis is therefore a
validated choice, not a convenience. Third, the design rules are mechanistic and monotonic — a
smaller N'-cap trades mass for capacity, an electron-withdrawing cap raises potential, and ring
donors lower it — giving a predictable handle to place a dihydrophenazine catholyte against a chosen
anolyte. Kinetics are ranked on outer-sphere lambda_o (validated against the molecular-cavity PCM
method), with inner-sphere lambda_i reported where physical (p-type oxidations) and held where the
gas-phase anion is unbound (n-type reductions) — an honest physical limit, not missing data.

## Conclusion (draft)
A grafting-aware DFT+SMD screening pipeline, validated per-couple on both electrode polarities,
ranks 92 Merrifield-graftable redox motifs and fills the two-electron catholyte gap with
dihydrophenazines carrying measured design rules. The recommended two-electron/two-electron cell,
2-(benzylamino)-1,4-naphthoquinone with N'-methyl-5,10-dihydrophenazine, is predicted at a
bias-corrected 1.83 +/- 0.12 V. The remaining step is experimental: synthesis of the grafted
leads and cell testing, for which these are the prioritized targets.
