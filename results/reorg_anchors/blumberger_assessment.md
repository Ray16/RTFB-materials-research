# Blumberger-group & related condensed-phase ET reorganization-energy anchors

**Purpose.** Build a high-trust validation anchor set for reorganization energy (λ) for the
redox-flow-battery redoxmer project, focused on the Blumberger group (Jochen Blumberger, UCL)
and closely related careful QM/MD condensed-phase electron-transfer (ET) studies of **organic
and molecular** systems.

**Companion file:** `blumberger_anchors.csv` (same directory).

## Retrieval note (important for trust)
Web *summarization* tools (WebSearch / WebFetch content step) were down for this session with a
backend fast-model error (`claude-haiku-4.5` unavailable). All citations and DOIs below were
therefore retrieved and verified through **primary bibliographic REST APIs** — OpenAlex,
Crossref, Europe PMC — and the open-access full text of the Blumberger 2017 review and the
Giannini 2019 paper was downloaded and parsed locally. **DOIs are verified** (resolve 200 on
Crossref). **Numeric λ values marked "lit-standard", "approx", "representative" or
"needs-verification" were taken from field knowledge / the primary compilation and were NOT
re-read digit-by-digit from the source PDF this session** — confirm exact digits before quoting
in a paper. Nothing here is fabricated; uncertain items are flagged as such.

## How our pipeline computes λ (for the reproducibility mapping)
- **Inner-sphere** `src/redox/reorg.py`: gas-phase **4-point Nelsen** (adiabatic potential
  energy surface) scheme, `λ_i = [E_O(q_R)−E_O(q_O)] + [E_R(q_O)−E_R(q_R)]`, using gas single
  points at **wB97M-V/def2-TZVP** on the SMD(MeCN)-optimized geometries. This is *exactly* the
  method Nelsen introduced (JACS 1987, 10.1021/ja00237a007) and the same 4-point definition
  used for organic-semiconductor internal λ throughout the literature.
- **Outer-sphere** `src/redox/solvated_reorg.py`: **Marcus/Born dielectric-continuum** term
  (Pekar factor `1/ε_op − 1/ε_r`, MeCN ε_r=37.5, ε_op=1.806) with a SASA-based effective radius;
  `λ_total = λ_i + λ_o`. This is a **liquid dielectric-continuum** model — not explicit-solvent
  MD and not a polarizable-crystal model.
- We do **not** run explicit-solvent / QM-MM MD energy-gap sampling.

## Trust tiers of the anchors

### Tier 1 — Inner-sphere (intramolecular) λ of organic molecules  →  DIRECTLY reproducible (a)
The most trustworthy *and* most reproducible anchors are the **internal (inner-sphere)
reorganization energies of the oligoacenes** (naphthalene, anthracene, tetracene, pentacene,
rubrene). These are computed by the **same gas-phase 4-point scheme our pipeline uses**, are
reported to ~10 meV consistency across many groups, and are the ones the Blumberger group adopts
for organic-semiconductor charge-transport MD.

- Primary numeric compilation: **Coropceanu, Cornil, da Silva Filho, Olivier, Silbey, Brédas,
  *Chem. Rev.* 2007, 107, 926 (10.1021/cr050140x)** — hole internal λ ≈ 0.18 (naphthalene),
  0.14 (anthracene), 0.11 (tetracene), 0.10 (pentacene), 0.15 (rubrene) eV, at ~B3LYP/6-31G(d,p)
  gas-phase 4-point. Electron internal λ are somewhat larger and less uniformly reported —
  treat electron values as needs-verification.
- Scale independently confirmed by the Blumberger methods review **Oberhofer, Reuter, Blumberger,
  *Chem. Rev.* 2017, 117, 10319 (10.1021/acs.chemrev.7b00086)**, which states organic-
  semiconductor λ ≈ **0.1 eV** (vs ≈ 1 eV for biological ET), and by **Giannini et al.,
  *Nat. Commun.* 2019, 10, 3843 (10.1038/s41467-019-11775-9)**, which states OSC λ is
  "**0.2 eV or less**" with a typical value of **150 meV**.

**How to reproduce:** run our inner-sphere 4-point on each neutral↔cation couple. To get a
*like-for-like* number set the level of theory to the literature level — **B3LYP/6-31G(d,p),
gas phase** — rather than our default wB97M-V/def2-TZVP (which will land close, within
~0.02–0.05 eV, but not identical). This is a genuine, cheap, high-value validation: these are
small rigid neutral molecules with no charge/spin ambiguity and no stereocenters.

**Bonus family (also Tier 1, method-exact):** **Nelsen's organic amine self-exchange couples**
(tetraalkylhydrazines, N,N,N′,N′-tetramethyl-*p*-phenylenediamine, etc.; Nelsen JACS 1987,
10.1021/ja00237a007) are the *origin* of the 4-point method — reproducing them is a direct test
of our own scheme. Specific per-molecule λ were not retrieved this session; worth pulling from
Nelsen's papers as an additional inner-sphere anchor set.

### Tier 2 — Outer-sphere / external λ of organic crystals  →  NOT directly reproducible (c)
- **McMahon, Troisi, *J. Phys. Chem. Lett.* 2010, 1, 941 (10.1021/jz1001049)** — external
  (outer-sphere) reorganization energy for hole transport in naphthalene→rubrene from a
  **polarizable force field for the molecular crystal**. Values are small (tens of meV; exact
  numbers needs-verification).
- Our outer-sphere term is a **dielectric-continuum liquid** model, physically different from a
  polarizable *solid*. So these crystal external-λ values are **not** reproducible by our
  `solvated_reorg.py` add-on. (Our continuum term is appropriate for *solution* outer-sphere λ,
  which is a different quantity.)

### Tier 3 — Total λ from explicit QM/MD (aqueous ions, proteins)  →  NOT reproducible (c)
These are the "highest-trust" values in the classical sense (full statistical-mechanical
sampling of the Marcus energy gap), but they are **total** λ dominated by explicit-solvent /
protein reorganization and require MD we do not run:
- **Aqueous transition-metal self-exchange** (Ru(H₂O)₆²⁺ᐟ³⁺, Fe²⁺ᐟ³⁺, Ru(bpy)₃²⁺ᐟ³⁺):
  total λ ≈ 1.5–2.2 eV (representative; per-system verification needed). Compiled/derived in
  **Blumberger, *Chem. Rev.* 2015, 115, 11191 (10.1021/acs.chemrev.5b00298)**; Ru(bpy)₃ via
  Seidel/Winter/Blumberger DFT-MD + liquid-microjet photoemission.
- **Protein ET**: cytochrome c (**Jiang, Futera, Blumberger, *J. Phys. Chem. B* 2019, 123, 7588,
  10.1021/acs.jpcb.9b05253**) and multiheme cytochromes (**Breuer, Rosso, Blumberger, Butt,
  *J. R. Soc. Interface* 2014, 12, 20141117, 10.1098/rsif.2014.1117**): total λ ≈ 0.6–1.5 eV
  (needs-verification; ergodicity/sampling-dependent).
- Not reproducible by our 4-point or continuum add-on. Useful only as *context* for the
  magnitude of outer-sphere contributions, not as a numeric target for us.

## Reproducibility verdict (bottom line)
- **What we can validate against directly (a):** the **inner-sphere (intramolecular) λ of neutral
  organic molecules** — the oligoacenes are the ready-made anchor set, and Nelsen's amine
  self-exchange couples are a method-exact bonus. Do it by setting our inner-sphere 4-point to
  **B3LYP/6-31G(d,p) gas phase** to match the literature level; expect agreement to ~0.01–0.05 eV.
  This is the single most defensible thing to do with our current code.
- **What needs an outer-sphere addition (b):** *solution-phase* total λ of a small ion where a
  dielectric-continuum outer-sphere term is defensible. None of the Blumberger organic anchors
  are cleanly in this bucket (their outer-sphere is crystal-polarizable or explicit-MD), so (b)
  is largely **not applicable** to this particular anchor list — flagged (c) instead where the
  environment is a crystal or protein.
- **What we cannot reproduce (c):** all **crystal external λ** and all **explicit-solvent/protein
  total λ** (aqueous metals, cytochromes) — these need polarizable-crystal or full QM/MM MD
  energy-gap sampling.

## Recommended level-of-theory settings to match anchors
- **Inner-sphere oligoacene benchmark:** B3LYP/6-31G(d,p), gas phase, 4-point Nelsen. (This is
  the historical level; a modern range-separated hybrid is arguably "better" physics but will
  *not* reproduce the tabulated numbers — for validation, match the source functional/basis.)
- Keep our production default (wB97M-V/def2-TZVP) for the project's own screening numbers, and
  report the B3LYP/6-31G(d,p) run only as the validation cross-check, so the anchor comparison is
  level-matched. (Analogous to the existing D3TaLES cross-check, which matches their IP-tuned
  LC-ωHPBE/def2-SVP level rather than ours.)

## Gotchas
- **Level-of-theory sensitivity:** internal λ shifts by ~0.02–0.05 eV between functionals/bases;
  a "disagreement" with a literature number is usually a level mismatch, not a bug. Always
  level-match before claiming (dis)agreement.
- **Inner vs total confusion:** literature "reorganization energy" for OSCs usually means the
  *internal* (gas-phase 4-point) part; for aqueous/biological ET it usually means the *total*
  (MD, outer-sphere-dominated). Do not compare our λ_i to a total-λ number.
- **Crystal ≠ solution outer-sphere:** our Born/Marcus continuum is a liquid-dielectric model;
  it is not the right physics for a molecular-crystal external λ (McMahon–Troisi). Don't use it
  to "reproduce" crystal external λ.
- **Diffuse functions for anions:** electron (anion) internal λ and any anion single point are
  sensitive to diffuse functions and to whether the gas-phase anion is *bound* — our pipeline
  already flags `anion_unbound` (gas HOMO > 0); for anchor electron-λ use a basis with diffuse
  functions (e.g. 6-31+G(d,p)) or score in SMD, consistent with that guard in `reorg.py`.
- **Numeric values here are pending exact-digit verification** (web summarizer was down); the
  DOIs are verified, so re-pull the exact tables from the primary sources when web tools recover.

## Verified citations (DOIs resolve on Crossref)
1. Coropceanu, Cornil, da Silva Filho, Olivier, Silbey, Brédas. "Charge Transport in Organic
   Semiconductors." *Chem. Rev.* 2007, 107, 926. doi:10.1021/cr050140x
2. Oberhofer, Reuter, Blumberger. "Charge Transport in Molecular Materials: An Assessment of
   Computational Methods." *Chem. Rev.* 2017, 117, 10319. doi:10.1021/acs.chemrev.7b00086
3. Blumberger. "Recent Advances in the Theory and Molecular Simulation of Biological Electron
   Transfer Reactions." *Chem. Rev.* 2015, 115, 11191. doi:10.1021/acs.chemrev.5b00298
4. Giannini, Carof, Ellis, Yang, Ziogos, Ghosh, Kubas, Blumberger. "Quantum localization and
   delocalization of charge carriers in organic semiconducting crystals." *Nat. Commun.* 2019,
   10, 3843. doi:10.1038/s41467-019-11775-9
5. McMahon, Troisi. "Evaluation of the External Reorganization Energy of Polyacenes." *J. Phys.
   Chem. Lett.* 2010, 1, 941. doi:10.1021/jz1001049
6. Jiang, Futera, Blumberger. "Ergodicity-Breaking in Thermal Biological Electron Transfer?
   Cytochrome c Revisited." *J. Phys. Chem. B* 2019, 123, 7588. doi:10.1021/acs.jpcb.9b05253
7. Breuer, Rosso, Blumberger, Butt. "Multi-haem cytochromes in Shewanella oneidensis MR-1:
   structures, functions and opportunities." *J. R. Soc. Interface* 2014, 12, 20141117.
   doi:10.1098/rsif.2014.1117
8. Nelsen, Blackstock, Kim. "Estimation of inner shell Marcus terms for amino nitrogen compounds
   by molecular orbital calculations." *J. Am. Chem. Soc.* 1987, 109, 677. doi:10.1021/ja00237a007
   (origin of the 4-point method our pipeline implements)
9. Kubas, Hoffmann, Heck, Oberhofer, Elstner, Blumberger. "Electronic couplings for molecular
   charge transfer: Benchmarking CDFT, FODFT, and FODFTB..." *J. Chem. Phys.* 2014, 140, 104105
   (HAB11; couplings, method context — not λ). erratum doi:10.1063/1.4916382
