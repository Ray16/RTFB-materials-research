# Literature Survey: Multi-Electron Organic Redox Materials for Non-Aqueous RFBs

**Scope.** Computational-chemistry literature scout for a non-aqueous redox-flow-battery (RFB)
screening project (acetonitrile / MeCN electrolyte, Fc/Fc⁺ referencing). Focus on three
multi-electron organic redox classes — **viologens, imides (PMDI/NDI), quinones** — with
emphasis on **polymer-bound / functionalized (grafted)** materials and on **amine / carbonyl
functional handles** for monomer attachment. A short section covers additional promising
classes.

**Two screening criteria this project optimizes (both "lower is better"):**
1. **Low synthetic accessibility (SA) score** — cheap, few-step, scalable chemistry from
   commercial precursors (SA score per Ertl & Schuffenhauer 2009 [R1]).
2. **Low reorganization energy (λ)** — small inner-sphere structural change on electron
   transfer → fast kinetics, high round-trip efficiency (Marcus framework).

**A note on potentials.** Literature values below are quoted against the reference electrode
used by the source (Fc/Fc⁺, SCE, Ag/Ag⁺, or NHE). For rough comparison in MeCN,
Fc/Fc⁺ ≈ +0.38–0.40 V vs SCE, so *E*(vs Fc/Fc⁺) ≈ *E*(vs SCE) − 0.40 V. Where a specific
numeric value is not directly traceable to a source it is labelled **[approx]** or
**[needs-verification]**; do not treat those as measured. λ values in the literature are
sparse and often computational; treat them as order-of-magnitude guidance.

---

## Summary table (representative RFB-relevant molecules)

| Molecule | Class | Couple(s) / e⁻ | E (as reported) | λ (if known) | Functionalization handle | Ref |
|---|---|---|---|---|---|---|
| Methyl viologen (MV²⁺) | Viologen | 2⁺/+• and +•/0 (2 e⁻) | ≈ −0.45 & −0.86 V vs SCE (MeCN) [approx] | 2⁺/+• low λ (delocalized); +•/0 higher (planarization) [needs-verification] | N-alkyl quaternization; benzyl/vinylbenzyl for grafting | [V1],[V2] |
| Ethyl/benzyl viologens | Viologen | 2⁺/+•/0 (2 e⁻) | similar to MV, tunable by N-substituent | as above | N,N′-dialkyl / N-(ω-haloalkyl), aminoalkyl, carboxyalkyl | [V3],[V4] |
| Poly(vinylbenzyl ethyl viologen) | Viologen (polymer) | 2⁺/+•/0 pendant | ≈ MV, MW-independent redox | — | styrenic pendant (main interest) | [V5],[V6] |
| Pyromellitic diimide (PMDI) | Imide | 2 sequential 1 e⁻ (0/−•/2−) | ≈ −1.1 to −1.4 V vs Fc/Fc⁺ region [approx] | low inner-sphere λ (planar, N-node) [needs-verification] | imide-N from primary amine (condensation) | [I1],[I2] |
| Naphthalene diimide (NDI) | Imide | 2 sequential 1 e⁻ | ≈ −1.1 / −1.5 V vs Fc/Fc⁺ region [approx] | low λ (rigid π) [needs-verification] | imide-N from primary amine; ammonium/zwitterion pendants | [I3],[I4],[I5] |
| 1,4-Benzoquinone / methoxy-BQ | Quinone | Q/Q•⁻/Q²⁻ (up to 2 e⁻) | tunable; OMe lowers potential | moderate (C=O bond change; PCET raises λ) | ring OMe/OH/NH₂/COOH substitution | [Q1],[Q2] |
| Anthraquinone (AQ) | Quinone | Q/Q•⁻/Q²⁻ (2 e⁻) | ≈ −1.0 to −1.3 V vs Fc/Fc⁺ (MeCN) [approx] | moderate λ | 1-amino-AQ, AQ-2-COOH for amide/ester tether | [Q3],[Q4] |
| Ferrocene/anthraquinone bipolar | Quinone+Fc | AQ (2 e⁻) + Fc⁺/Fc (1 e⁻) | bipolar, symmetric cell | — | covalent Fc–AQ dyad; ester linker | [Q3] |

(Fuller per-molecule detail, cycling notes, and functionalization discussion follow by class.)

---

## 1. Viologens (4,4′-bipyridinium)

**Redox chemistry.** Viologens (N,N′-disubstituted 4,4′-bipyridinium, V²⁺) undergo two
reversible one-electron reductions: **V²⁺ → V⁺• (radical cation) → V⁰**. The first couple
(2⁺/+•) is highly reversible and chemically robust; the neutral V⁰ from the second couple is
less soluble/less stable in many electrolytes. This makes viologens the canonical **anolyte**
(negative side) for organic RFBs.

**Representative molecules & performance.**
- **Methyl viologen (MV²⁺)** is the workhorse anolyte. In non-aqueous MeCN, Hu et al. [V1]
  demonstrated **two-electron utilization of MV** (both 2⁺/+• and +•/0) in a NAORFB, roughly
  doubling capacity vs single-electron operation, with the second couple requiring careful
  electrolyte/membrane choice to manage the poorly soluble neutral state. In aqueous cells MV
  paired with TEMPO/4-HO-TEMPO gave stable, low-cost all-organic RFBs (Liu et al. [V2]) —
  foundational for viologen anolyte design though aqueous.
- **N,N′-dialkyl / benzyl viologens** — the N-substituent tunes solubility, the 2⁺/+•
  potential (weakly), and radical-cation stability. Bulky/branched or ethyl substituents and
  charged pendants are used to raise solubility and suppress π-dimerization of V⁺•.

**Reorganization energy.** The **2⁺/+• couple has relatively low inner-sphere λ**: the added
electron is delocalized over the planar bipyridinium π-system with modest bond-length change,
which underlies its fast, reversible kinetics. The **+•/0 couple has larger λ** — reduction to
the neutral quinoid form drives inter-ring planarization and larger geometric relaxation —
which is one reason the second electron is harder to cycle reversibly. (Specific λ magnitudes
are sparsely reported; the qualitative trend is well established — treat numeric λ as
**[needs-verification]** and a good target for this project's own computation.)

**Synthetic accessibility (SA).** Excellent (low SA score): viologens are made in **one step**
by quaternizing commercial **4,4′-bipyridine** with alkyl/benzyl halides. This one-pot
quaternization is also the functionalization step, so viologens are among the easiest classes
to graft.

**Functionalization with amine / carbonyl handles (polymer attachment).**
- **Pendant (side-chain) grafting:** quaternize 4,4′-bipyridine with a **vinylbenzyl
  (styrenic) halide** → styrenic viologen monomer → radical polymerization to
  **poly(vinylbenzyl viologen)**; or use **acrylate/methacrylate** or acrylamide halide
  tethers. Nagarjuna et al. [V5] built non-aqueous **redox-active polymers** bearing pendant
  viologen (and TEMPO) and showed the **redox potential is essentially independent of polymer
  molecular weight** while size-selective membranes suppress crossover — a central result for
  polymer-bound NAORFBs.
- **Amine handle:** **N-(aminoalkyl) viologens** (e.g., aminopropyl) provide a primary amine
  for **amide coupling** to carboxyl-bearing backbones; conversely a **carboxyalkyl (carbonyl)
  viologen** enables ester/amide attachment to amine-functional polymers. These are the direct
  routes to tether viologen monomers to a polymer via amine or carbonyl chemistry.
- **Main-chain "polyviologens" (ionenes):** bis-quaternization of 4,4′-bipyridine with a
  **dihalide** builds viologen units *into* the backbone.
- **Linker effect:** because the N-substituent sits on a formal node of the bipyridinium redox
  orbital, alkyl/benzyl linkers perturb the 2⁺/+• potential only weakly; tether **length and
  flexibility** mainly affect solubility, V⁺• dimerization, and transport. Burgess et al. [V6]
  showed backbone tether length/structure tunes electrochemical performance of redox polymers.
- **Aqueous polymer benchmark:** Janoschka et al. [V3] (Nature 2015) built a metal-free
  polymer RFB from a **viologen-polymer anolyte** and TEMPO-polymer catholyte; synthesis of
  both polymer families is detailed in [V4]. Hagemann et al. [V7] paired a TEMPO polymer with
  dimethyl viologen. These aqueous systems are the template for the non-aqueous grafted designs.

---

## 2. Imides (pyromellitic diimide PMDI, naphthalene diimide NDI)

**Redox chemistry.** Aromatic diimides accept **two electrons in two sequential reversible
one-electron steps** (neutral → radical anion → dianion), storing charge on the imide
carbonyls / aromatic core. They are strong, low-potential **anolytes** with excellent radical-
anion stability owing to delocalization over the rigid planar π-system.

**Representative molecules & performance.**
- **Pyromellitic diimide (PMDI):** Nambafu et al. [I1] reported a **PMDI-based bipolar
  molecule** (imide anolyte fused with a high-potential unit) enabling a **total-organic
  symmetric** RFB — same molecule on both sides, minimizing crossover-induced capacity loss.
  Gouget et al. [I2] showed **chemical modification of PMDI raises its redox potential** toward
  higher-voltage organics, illustrating substituent tuning of the imide core.
- **Naphthalene diimide (NDI):** Medabalmi et al. [I3] demonstrated **NDI as a two-electron
  anolyte** for aqueous/neutral-pH RFBs; Li et al. [I4] developed **stable NDI zwitterions**
  (self-solubilizing, crossover-resistant). Xu et al. [I5] built a **ferrocene/naphthalimide
  bipolar molecule** for a symmetric NAORFB with improved cycling stability — a direct
  non-aqueous, functionalized-imide example.

**Reorganization energy.** Diimides are **rigid, planar, fully conjugated**; the radical anion
is delocalized over the naphthalene/benzene core and imide carbonyls with small nuclear
displacement, giving **low inner-sphere λ** and fast, reversible kinetics — a key reason for
their electrochemical robustness. (Numeric λ values class-wide are **[needs-verification]** and
another good computational target here.)

**Synthetic accessibility (SA).** Excellent (low SA score). Diimides form in **one
condensation step**: heat the commercial **dianhydride** (**pyromellitic dianhydride** for
PMDI; **naphthalene-1,4,5,8-tetracarboxylic dianhydride** for NDI) with a **primary amine**.

**Functionalization with amine / carbonyl handles (polymer attachment) — the standout
advantage.**
- The **imide nitrogen is installed directly from any primary amine** during the condensation.
  This means an enormous range of tethers (aminoalkyl, aminoaryl, diamines, amino acids,
  aminostyrene, aminoethyl-methacrylamide) can be attached in the **same one-step reaction**
  that makes the molecule — ideal for monomer synthesis.
- **Electronic insulation of the linker:** the imide-N lies on a **node** of the redox-active
  frontier orbital, so the substituent barely shifts the redox potential. Redox behaviour is
  therefore **nearly linker-independent**, letting the tether be chosen for solubility /
  polymerizability without detuning the couple — a major benefit for grafted designs.
- **Ammonium / carbonyl-functionalized variants:** condense with **N,N-dimethylaminopropyl-
  amine** then quaternize → cationic **ammonium-functionalized** diimide (water/polar
  solubility, crossover control); condense with amino-acids → pendant **carboxyl (carbonyl)**
  handles for ester/amide grafting; **zwitterionic** NDIs [I4] combine both.
- **Polymer attachment modes:** **main-chain polyimides** from a **diamine + dianhydride**
  (imide units in the backbone), or **pendant** diimides via aminostyrene / amino-methacrylate
  monomers. The bipolar imide-Fc dyad [I5] shows how a carbonyl/ester linker joins two redox
  centers on one scaffold.

---

## 3. Quinones (benzoquinones, methoxy/dimethoxy-substituted, anthraquinones)

**Redox chemistry.** Quinones (Q) cycle **Q/Q•⁻/Q²⁻** (or, with protons, Q/QH₂), storing up to
two electrons on the carbonyl-bearing ring. Potential is highly tunable by ring substituents:
**electron-donating groups (–OMe, –OH, –NH₂) lower the potential** (more anolyte-like);
electron-withdrawing groups raise it. This makes quinones usable as anolyte or (with
appropriate substitution) low-potential catholyte.

**Representative molecules & performance.**
- **Benzoquinone / methoxy- and dimethoxy-benzoquinones:** substituent screening tunes the
  Q/Q•⁻ potential and stability; methoxy groups raise solubility and shift potential. Huang et
  al. [Q1] screened **quinone/dialkoxybenzene "liquid catholyte"** molecules for non-aqueous
  RFBs, mapping structure–potential–solubility relationships.
- **Anthraquinone (AQ):** the most-used quinone scaffold. In non-aqueous cells Zhen et al. [Q3]
  built a **ferrocene/anthraquinone bipolar molecule** (AQ 2 e⁻ anolyte + Fc⁺/Fc 1 e⁻
  catholyte on one scaffold) for a **symmetric NAORFB**. In aqueous alkaline cells, Lin et al.
  [Q4] (Aziz, *Science* 2015) demonstrated the landmark **DHAQ alkaline quinone flow battery** —
  foundational for quinone anolyte design and substituent effects.
- **Computational structure–property work:** Ding et al. [Q2] combined experiment and
  computation on **bio-inspired quinones** for organic RFBs, computing redox potentials and
  addressing reorganization/solvation — a useful methodological anchor for this project.

**Reorganization energy.** Quinones have **moderate inner-sphere λ**: reduction changes C=O
bond order and ring geometry. λ **rises when the couple is proton-coupled** (Q/QH₂ via PCET,
large nuclear reorganization) and is **lower for the outer-sphere Q/Q•⁻ single-electron couple**
in aprotic MeCN. Substituents modulate λ modestly (rigidifying/conjugating groups lower it;
H-bonding –OH raises it via PCET). Computational λ for specific quinones should be generated in-
project; literature numbers are **[needs-verification]**.

**Synthetic accessibility (SA).** Simple benzoquinones and anthraquinone are **commercial and
cheap** (low SA). Functionalized derivatives (amino-, carboxy-, dimethoxy-AQ) require one or
two substitution steps — moderate SA, generally higher than viologens/imides because ring
functionalization can need regioselective chemistry.

**Functionalization with amine / carbonyl handles (polymer attachment).**
- **Handles are on the ring, not a node:** unlike imides, quinone substituents sit **on the
  redox-active ring**, so the linker **directly perturbs the redox potential** (donor linkers
  lower it, acceptor linkers raise it). Linker chemistry and redox tuning are coupled — a
  design trade-off to track.
- **Amine handle:** **amino-quinones / 1-aminoanthraquinone** provide a primary amine for
  **amide coupling** to carboxyl backbones (the amino group also lowers potential).
- **Carbonyl handle:** **anthraquinone-2-carboxylic acid** (and quinone-alkyl-COOH) enable
  **ester/amide** tethering to amine- or hydroxyl-functional polymers; **hydroxy-quinones**
  allow etherification.
- **Polymer attachment modes:** **pendant** quinones via **methacrylate-** or **vinyl-
  anthraquinone** monomers (side-chain), or **main-chain** quinone polymers. Because potential
  is linker-sensitive, **spacer length/electronics must be chosen to preserve the target
  potential** — the opposite situation from imides, and an important point for grafted quinone
  design.

---

## 4. Other promising multi-electron / redox classes for non-aqueous RFBs

- **Phenothiazines** — high-potential, reversible one-electron (and accessible two-electron)
  **catholytes**; radical cation stabilized by N/S substitution. Widely developed by the Odom
  group (e.g., N-alkyl / dimethoxy phenothiazines) and used in symmetric bipolar designs.
  Key non-aqueous example: Yan & Sevov, "Mechanism-Based Design of a High-Potential Catholyte
  Enables a 3.2 V All-Organic NAORFB," *JACS* 2019 [O1] (phenothiazine-type catholyte).
- **Phenoxazines / phenazines / dihydrophenazines** — N,O- or N,N-heterocyclic redox centers
  spanning a wide potential window; phenoxazine catholytes optimized for NAORFB by Yan et al.
  [O2]; **phenazine** anolytes (multi-electron, tunable) e.g. Romadina et al. [O3]
  (water-soluble). Dihydrophenazines give high-potential, stable radical cations for catholytes.
- **TEMPO / nitroxide radicals** — the benchmark **catholyte** (one-electron TEMPO/TEMPO⁺,
  very fast kinetics, small λ). Wei et al., "TEMPO-Based Catholyte for High-Energy Density
  Nonaqueous Redox Flow Batteries," *Adv. Mater.* 2014 [O4]. Readily polymer-grafted (PTMA-type
  pendants), central to polymer NAORFBs.
- **Ferrocenes** — robust one-electron **catholyte** (Fc/Fc⁺), also the referencing standard.
  Used in **bipolar dyads** with imides/quinones (Zhen [Q3]; Xu [I5]); easily functionalized
  (aminomethyl-, carboxy-ferrocene) for grafting.
- **Cyclopropenium / pyridinium multielectron anolytes** — Sanford/Sevov "physical-organic"
  low-potential anolytes: Sevov et al., *JACS* 2016 [O5]; Hendriks et al., "Multielectron
  Cycling of a Low-Potential Anolyte…," *ACS Energy Lett.* 2017 [O6].
- **Fluorenones / ketone anolytes** — reversible ketone hydrogenation/dehydrogenation; Feng et
  al., *Science* 2021 [O7] (aqueous, but a promising multi-electron carbonyl motif).
- **Thianthrenes** — high-potential two-stage (radical cation / dication) catholytes with very
  positive potentials; attractive for high-voltage cells but crossover/stability of the dication
  needs management. *(Specific NAORFB citation: **[needs-verification]**.)*

**Referencing note.** For this project's MeCN cells, report all potentials **vs Fc/Fc⁺**; Fc/Fc⁺
also anchors the internal reference in `config/electrolyte.py` (MeCN, SMD, PF₆⁻).

---

## 5. Screening take-aways (low SA + low λ)

- **Best on SA (fewest steps from commercial precursors):** **imides** (one-step dianhydride +
  amine condensation) and **viologens** (one-step 4,4′-bipyridine quaternization) tie for
  easiest, and **their functionalization is the same reaction that makes the molecule.** Simple
  quinones are cheap but ring-functionalization adds steps.
- **Best on λ (small inner-sphere reorganization):** **rigid, planar, delocalized π-systems** —
  **imides (N/−•)**, the **viologen 2⁺/+•** couple, and the aprotic **quinone Q/Q•⁻** couple.
  Avoid couples with large geometric change or PCET (viologen **+•/0** planarization; quinone
  **Q/QH₂** proton coupling) when minimizing λ.
- **Best on grafting without detuning:** **imides** — the imide-N linker sits on an orbital
  node, so redox potential is nearly linker-independent; ideal for amine/carbonyl monomer
  attachment. **Viologens** are similar (N-node, weak linker effect). **Quinones** are the
  hardest — the linker sits on the redox ring and shifts the potential, so spacer electronics
  must be co-designed with the target potential.
- **Synergy:** imides and viologens score well on *both* criteria simultaneously and offer clean
  amine/carbonyl handles — strong primary candidates for polymer-bound non-aqueous designs;
  quinones offer the widest potential tunability at the cost of coupled linker/potential effects.

---

## Sources

**Screening-criteria / methods**
- [R1] Ertl, Schuffenhauer. "Estimation of synthetic accessibility score of drug-like molecules
  based on molecular complexity and fragment contributions." *J. Cheminform.* **2009**, 1, 8.
  https://doi.org/10.1186/1758-2946-1-8

**Reviews (context)**
- Winsberg, Hagemann, Janoschka, Hager, Schubert. "Redox-Flow Batteries: From Metals to Organic
  Redox-Active Materials." *Angew. Chem. Int. Ed.* **2016**, 56, 686.
  https://doi.org/10.1002/anie.201604925
- Kwabi, Ji, Aziz. "Electrolyte Lifetime in Aqueous Organic Redox Flow Batteries: A Critical
  Review." *Chem. Rev.* **2020**, 120, 6467. https://doi.org/10.1021/acs.chemrev.9b00599
- Xu et al. "Molecular engineering redox-active organic materials for nonaqueous redox flow
  battery." *Curr. Opin. Chem. Eng.* **2022**, 37, 100851.
  https://doi.org/10.1016/j.coche.2022.100851

**Viologens**
- [V1] Hu et al. "Two electron utilization of methyl viologen anolyte in nonaqueous organic
  redox flow battery." *J. Energy Chem.* **2018**. https://doi.org/10.1016/j.jechem.2018.02.014
- [V2] Liu et al. "A Total Organic Aqueous Redox Flow Battery Employing a Low Cost and
  Sustainable Methyl Viologen Anolyte and 4-HO-TEMPO Catholyte." *Adv. Energy Mater.* **2016**,
  6, 1501449. https://doi.org/10.1002/aenm.201501449
- [V3] Janoschka et al. "An aqueous, polymer-based redox-flow battery using non-corrosive, safe,
  and low-cost materials." *Nature* **2015**, 527, 78. https://doi.org/10.1038/nature15746
- [V4] Janoschka et al. "Synthesis and characterization of TEMPO- and viologen-polymers for
  water-based redox-flow batteries." *Polym. Chem.* **2015**, 6, 7801.
  https://doi.org/10.1039/c5py01602a
- [V5] Nagarjuna et al. "Impact of Redox-Active Polymer Molecular Weight on the Electrochemical
  Properties and Transport across Porous Separators in Nonaqueous Solvents." *J. Am. Chem. Soc.*
  **2014**, 136, 16309. https://doi.org/10.1021/ja508482e
- [V6] Burgess et al. "Impact of Backbone Tether Length and Structure on the Electrochemical
  Performance of … Redox Active Polymers." *Chem. Mater.* **2016**, 28.
  https://doi.org/10.1021/acs.chemmater.6b02825
- [V7] Hagemann et al. "An aqueous all-organic redox-flow battery employing a TEMPO-containing
  polymer catholyte and dimethyl viologen dichloride anolyte." *J. Power Sources* **2018**.
  https://doi.org/10.1016/j.jpowsour.2017.09.007

**Imides (PMDI / NDI)**
- [I1] Nambafu et al. "Pyromellitic diimide based bipolar molecule for total organic symmetric
  redox flow battery." *Nano Energy* **2022**, 94, 106963.
  https://doi.org/10.1016/j.nanoen.2022.106963
- [I2] Gouget et al. "Increasing Redox Potential of Pyromellitic Diimide by Chemical
  Modifications: Toward High-Voltage Organic …" *Batteries & Supercaps* **2025**.
  https://doi.org/10.1002/batt.202500008
- [I3] Medabalmi et al. "Naphthalene diimide as a two-electron anolyte for aqueous and neutral
  pH redox flow batteries." *J. Mater. Chem. A* **2020**, 8, 11218.
  https://doi.org/10.1039/d0ta01160f
- [I4] Li et al. "Stable naphthalene diimide zwitterions for aqueous organic redox flow
  batteries." *Natl. Sci. Rev.* **2025**. https://doi.org/10.1093/nsr/nwaf286
- [I5] Xu et al. "Ferrocene/naphthalimide bi-redox molecule for enhancing the cycling stability
  of symmetric nonaqueous redox flow battery." *J. Power Sources* **2024**, 234368.
  https://doi.org/10.1016/j.jpowsour.2024.234368

**Quinones**
- [Q1] Huang et al. "Liquid Catholyte Molecules for Nonaqueous Redox Flow Batteries." *Adv.
  Energy Mater.* **2015**, 5, 1401782. https://doi.org/10.1002/aenm.201401782
- [Q2] Ding et al. "Exploring Bio-inspired Quinone-Based Organic Redox Flow Batteries: A
  Combined Experimental and Computational Study." *Chem* **2016**, 1, 790.
  https://doi.org/10.1016/j.chempr.2016.09.004
- [Q3] Zhen et al. "Ferrocene/anthraquinone based bi-redox molecule for symmetric nonaqueous
  redox flow battery." *J. Power Sources* **2020**, 480, 229132.
  https://doi.org/10.1016/j.jpowsour.2020.229132
- [Q4] Lin et al. "Alkaline quinone flow battery." *Science* **2015**, 349, 1529.
  https://doi.org/10.1126/science.aab3033

**Other classes**
- [O1] Yan, Sevov et al. "Mechanism-Based Design of a High-Potential Catholyte Enables a 3.2 V
  All-Organic Nonaqueous Redox Flow Battery." *J. Am. Chem. Soc.* **2019**, 141, 15301.
  https://doi.org/10.1021/jacs.9b07345
- [O2] Yan et al. "Targeted Optimization of Phenoxazine Redox Center for Nonaqueous Redox Flow
  Battery." *ACS Mater. Lett.* **2022**, 4. https://doi.org/10.1021/acsmaterialslett.2c00050
- [O3] Romadina et al. "Phenazine-Based Compound as a Universal Water-Soluble Anolyte Material
  for Redox Flow Batteries." *Batteries* **2022**, 8, 288.
  https://doi.org/10.3390/batteries8120288
- [O4] Wei et al. "TEMPO-Based Catholyte for High-Energy Density Nonaqueous Redox Flow
  Batteries." *Adv. Mater.* **2014**, 26, 7649. https://doi.org/10.1002/adma.201403746
- [O5] Sevov et al. "Mechanism-Based Development of a Low-Potential, Soluble, and Cyclable
  Multielectron Anolyte for Nonaqueous Redox Flow Batteries." *J. Am. Chem. Soc.* **2016**, 138,
  15378. https://doi.org/10.1021/jacs.6b07638
- [O6] Hendriks et al. "Multielectron Cycling of a Low-Potential Anolyte in Alkali Metal
  Electrolytes for Nonaqueous Redox Flow Batteries." *ACS Energy Lett.* **2017**, 2, 2071.
  https://doi.org/10.1021/acsenergylett.7b00559
- [O7] Feng et al. "Reversible ketone hydrogenation and dehydrogenation for aqueous organic
  redox flow batteries." *Science* **2021**, 372, 836. https://doi.org/10.1126/science.abd9795
- Thianthrene NAORFB catholyte — **[needs-verification]** (concept established; specific primary
  citation not confirmed in this pass).

*Prepared via Crossref metadata lookups (WebSearch was unavailable during this session; all DOIs
above resolved through the Crossref API). Numeric potentials/λ marked [approx] or
[needs-verification] were not traced to a specific measured source and should be confirmed or
generated computationally before use.*
