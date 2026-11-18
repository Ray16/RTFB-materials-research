# Introduction

Non-aqueous organic redox flow batteries (ORFBs) are an attractive route to
grid-scale electrochemical storage because they decouple energy and power,
access cell voltages well beyond the aqueous stability window, and are built
from earth-abundant, structurally tunable organic redox couples rather than
scarce transition metals [Soloveichik 2015]. One of the more promising device
architectures immobilizes the redox-active species on a soluble or swellable
polymer, so that the active material is retained by a size-selective membrane
and crossover is suppressed [Janoschka 2015]. In practice these materials are
most readily made not by polymerizing a bespoke redox monomer but by grafting a
preformed redox motif onto a reactive, preformed backbone — most commonly a
chloromethylated (Merrifield-type) polystyrene, whose benzyl chloride handles
are alkylated by a nucleophilic site on the motif [Merrifield 1963]. This grafting
route decouples the chemistry of the backbone from the chemistry of the redox
center and makes a large library of candidate motifs synthetically plausible
from a single polymer platform. The design space it opens, however, is badly
unbalanced. The catalog of well-behaved organic anolytes — n-type couples based
on quinones, imides, and related electron-poor scaffolds — is deep and still
growing [Ding 2018], whereas robust, high-potential catholytes remain scarce
[Kwon 2019]. Because a flow cell needs both a catholyte and an anolyte, this
asymmetry, not the anolyte inventory, is what limits achievable cell voltage,
and it is precisely the kind of gap that a well-posed computational screen
should be able to survey and help fill.

Computation is in principle well suited to this task: redox potentials, and the
reorganization energetics that govern electron-transfer kinetics, are
accessible from solvated electronic-structure calculations, and density
functional theory (DFT) with an implicit solvent model has been used
extensively to rank organic redox couples [Marenich 2014; Er 2015]. In practice, though,
most published screens are difficult to trust or to reuse for the specific
problem of grafted motifs. Three shortcomings recur. First, the molecular model
is often a proxy for convenience rather than for the real material: a bare
parent heterocycle, an N–H or N–methyl analogue, or a soluble monomer stands in
for the grafted species, even though the attachment chemistry changes the
electronic environment at the redox center and therefore the potential. A capped
graft, a soluble parent, and an oligomer are genuinely different models and
should not be conflated. Second, reliability is rarely gated per candidate.
Spin contamination, unconverged or distorted geometries, missing charge states,
and inconsistent reference points can silently corrupt a reported potential,
and screens seldom document a per-state quality-control (QC) policy or state
which candidates were held back because they failed it. Third, and most
importantly, computed numbers are frequently reported without an external
calibration against measured potentials under matched conditions, so the reader
cannot tell whether a ranking reflects chemistry or method artifacts. Underlying
all three is a reproducibility problem: level of theory, solvent, standard
state, reference electrode, and software versions are often underspecified, and
the raw records needed to regenerate a figure are not released, so neither the
ranking nor its error bars can be independently reconstructed [Lejaeghere 2016]. A
computed descriptor, or a completed calculation, is not by itself evidence of
synthetic feasibility or of material performance; a screen is only as useful as
the discipline with which these distinctions are maintained.

This work addresses those gaps with a single, validated, end-to-end pipeline for
screening graftable redox motifs for non-aqueous ORFBs. Its design choices map
directly onto the shortcomings above. (i) Each motif is modeled as an *exact
capped graft*: the redox center carries a 4-methylbenzyl group as an explicit
Merrifield proxy, and the relevant oxidized or reduced charge states are built
and optimized individually, so the computed species is the grafted motif rather
than a stand-in. (ii) Every candidate passes explicit per-state QC —
spin/⟨S²⟩ checks, geometry integrity, full charge-state coverage, and
cross-point consistency — and candidates that fail are held rather than silently
dropped. (iii) All potentials are computed at ωB97M-V/def2-TZVPD with SMD(MeCN)
and reported against the ferrocene (Fc) reference in a TBAPF₆/MeCN electrolyte,
and the pipeline is calibrated against experiment, reproducing measured redox
potentials with a mean absolute error of 0.247 V and a systematic +0.19 V bias;
because the bias is systematic, relative rankings are more reliable than
absolute values. (iv) Candidates are ranked on a transparent, uncertainty-aware
(σ-aware) Pareto front over separate physical properties — potential and usable
capacity — with anolyte and catholyte pools treated separately and half-reaction
and reorganization quantities kept distinct; we deliberately avoid any
weighted-sum scalar score. Reorganization energetics are reported with
property-specific confidence: the outer-sphere term is a defensible continuum
estimate, whereas the inner-sphere term is currently treated as qualitative.
(v) The entire pipeline is byte-exact reproducible, distributed with portable
replays, a large contract-test suite, and a certificate-locked thermochemistry
kernel, so any reported number can be regenerated from its raw records.

We demonstrate the pipeline on a pool of 91 graftable candidates spanning
quinone, imide, pyridinic, phenazine, and phenothiazine families. The screen
confirms quantitatively that the accessible space is strongly anolyte-dominated
(85 anolyte versus 6 catholyte candidates), and the sparsity of the catholyte
Pareto front motivates a focused, matched-substituent study of grafted
5,10-dihydrophenazine (DMPZ) two-electron catholytes. From matched series
differing by a single substituent we extract an interpretable design rule —
electron-withdrawing N′-substituents raise the potential while ring-methylation
lowers it — that populates the previously empty catholyte front with
DFT-supported candidates near +0.3 to +0.6 V versus Fc. Pairing these catholytes
with pipeline anolytes suggests non-aqueous full cells near ~1.9 V; we stress
that these are computational proposals, not demonstrations of capacity, cycling,
crossover, or synthetic feasibility. The remainder of the paper is organized as
follows. We first describe the grafting abstraction and the nine-stage pipeline,
including its two-tier thermochemistry and its reliability-gating policy. We then
present the external validation against experiment that anchors the method's
credibility, apply the pipeline to the 91-candidate pool and its σ-aware Pareto
fronts, and finally present the DMPZ catholyte showcase and the proposed full
cells, before closing with the method's limitations and outlook.

## References (Introduction)

*DOI-verified against Crossref (2026-11-17). Author-year keys used in text; JCIM
numbering applied at submission.*

- **[Soloveichik 2015]** Soloveichik, G. L. Flow Batteries: Current Status and Trends. *Chem. Rev.* **2015**, *115* (20), 11533–11558. DOI: 10.1021/cr500720t.
- **[Janoschka 2015]** Janoschka, T.; Martin, N.; Martin, U.; Friebe, C.; Morgenstern, S.; Hiller, H.; Hager, M. D.; Schubert, U. S. An Aqueous, Polymer-Based Redox-Flow Battery Using Non-Corrosive, Safe, and Low-Cost Materials. *Nature* **2015**, *527* (7576), 78–81. DOI: 10.1038/nature15746.
- **[Merrifield 1963]** Merrifield, R. B. Solid Phase Peptide Synthesis. I. The Synthesis of a Tetrapeptide. *J. Am. Chem. Soc.* **1963**, *85* (14), 2149–2154. DOI: 10.1021/ja00897a025.
- **[Ding 2018]** Ding, Y.; Zhang, C.; Zhang, L.; Zhou, Y.; Yu, G. Molecular Engineering of Organic Electroactive Materials for Redox Flow Batteries. *Chem. Soc. Rev.* **2018**, *47* (1), 69–103. DOI: 10.1039/C7CS00569E.
- **[Kwon 2019]** Kwon, G.; Lee, K.; Lee, M. H.; Lee, B.; Lee, S.; Jung, S.-K.; Ku, K.; Kim, J.; Park, S. Y.; Kwon, J. E.; Kang, K. Bio-Inspired Molecular Redesign of a Multi-Redox Catholyte for High-Energy Non-Aqueous Organic Redox Flow Batteries. *Chem* **2019**, *5* (10), 2642–2656. DOI: 10.1016/j.chempr.2019.07.006.
- **[Marenich 2014]** Marenich, A. V.; Ho, J.; Coote, M. L.; Cramer, C. J.; Truhlar, D. G. Computational Electrochemistry: Prediction of Liquid-Phase Reduction Potentials. *Phys. Chem. Chem. Phys.* **2014**, *16* (29), 15068–15106. DOI: 10.1039/C4CP01572J.
- **[Er 2015]** Er, S.; Suh, C.; Marshak, M. P.; Aspuru-Guzik, A. Computational Design of Molecules for an All-Quinone Redox Flow Battery. *Chem. Sci.* **2015**, *6* (2), 885–893. DOI: 10.1039/C4SC03030C.
- **[Lejaeghere 2016]** Lejaeghere, K.; et al. Reproducibility in Density Functional Theory Calculations of Solids. *Science* **2016**, *351* (6280), aad3000. DOI: 10.1126/science.aad3000.
