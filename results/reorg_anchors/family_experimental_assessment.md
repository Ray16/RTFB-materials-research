# Experimental / literature reorganization-energy anchors for four redox families

Reproducibility assessment for the RFB redoxmer validation set. Companion data:
`results/reorg_anchors/family_experimental_anchors.csv`.

**Scope:** viologens, TEMPO/nitroxides, phenothiazines, quinones. Goal: anchor the
project's computed reorganization energy (λ) against experiment and against well-defined
computed references, and state exactly how each anchor maps onto our pipeline.

**Session limitation (read first).** The primary experimental papers for these families are
old and paywalled; the web-search/fetch backend was down this session (haiku model 400
errors). I reached Crossref/OpenAlex/Unpaywall/Semantic-Scholar **metadata** APIs by direct
HTTP, so every DOI below is verified to exist with the correct title/authors/year, **but the
numeric λ values inside those papers could not be extracted** (no OA full text; publisher
PDFs 403/JS-gated). Those rows have `lambda_eV` left blank and flagged "RETRIEVE". I did
**not** invent numbers. The quantitative backbone of this anchor set is therefore the
**computed** λ (D3TaLES + our own pipeline), which are real values I read out of the repo.

---

## 1. What our pipeline computes (the thing being validated)

- `src/redox/reorg.py` — **inner-sphere λ_i**, gas-phase, 4-point Nelsen scheme, at
  wB97M-V/def2-TZVP single points on SMD-optimized geometries. λ_i = geometric/vibrational
  distortion only. This is what most rows must be compared against.
- `src/redox/solvated_reorg.py` — **outer-sphere λ_o** via the Marcus/Born continuum
  (Pekar factor 1/ε_op − 1/ε_s), with an RDKit SASA hydrodynamic radius. `λ_total = λ_i + λ_o`.
- Electrolyte (`config/electrolyte.py` / `config/project.json`): **acetonitrile**, ε_s = 37.5,
  ε_op = 1.806 (Pekar factor = 0.527), Fc/Fc⁺ reference.

**Key fact for comparison:** experimental λ from a Marcus rate fit (self-exchange or
heterogeneous ET) or from an IV-CT band is almost always **TOTAL** (inner + outer). Our
default λ_i is **inner only**. To compare, add an outer-sphere term.

### Outer-sphere λ_o in MeCN (computed here, this session, CPU-only)

Two-sphere Marcus self-exchange λ_o at contact (= the repo's single-ion Born expression with
z=1); z=2 shown for the viologen dication step.

| Molecule | r_hyd (Å, SASA) | λ_o (1e) eV | λ_o (2e) eV |
|---|---|---|---|
| 1,4-benzoquinone | 4.01 | 0.95 | 3.79 |
| 1,4-naphthoquinone | 4.53 | 0.84 | 3.35 |
| 9,10-anthraquinone | 5.00 | 0.76 | 3.04 |
| duroquinone | 4.94 | 0.77 | 3.07 |
| TEMPO | 4.84 | 0.78 | 3.14 |
| phenothiazine | 5.01 | 0.76 | 3.03 |
| N-methylphenothiazine | 5.10 | 0.74 | 2.97 |
| methyl viologen (2+) | 5.30 | 0.72 | 2.87 |

**Caveats on λ_o.** (i) These are large (0.7–0.95 eV) and typically **dominate** total λ in
MeCN — so a 0.4 eV inner-sphere λ_i can correspond to a ~1.3 eV experimental total. (ii) The
radius is a structural SASA proxy, not a diffusion-measured r_hyd, so absolute λ_o is
approximate (±0.1–0.2 eV); trends are robust. (iii) **Homogeneous self-exchange vs
heterogeneous (electrode):** the electrode half-reaction gets an image-charge factor, so
λ_o,het ≈ ½ λ_o,homo. Match the geometry of the experiment you are comparing to. (iv) For the
**viologen 2+/+• step** do NOT use the z=2 column as "λ_o of the couple" — one electron is
transferred (Δz=1); the z² scaling applies to the Born solvation of a change in net charge,
not to the electrons transferred in the ET step. Use the 1e column for every couple here.

---

## 2. Per-family reproducibility verdict

### Quinones — STRONGEST experimental footing
- **Direct self-exchange (Marcus rate fits):** Williams 1969 (NMR line-broadening) gives
  **p-benzoquinone and duroquinone** molecule/radical-anion self-exchange rate constants —
  the two cleanest parent-quinone anchors. Grampp/Landgraf/Rasmussen 1999 gives **DDQ** in a
  MeCN solvent series with a solvent-dynamical (λ_o-separating) treatment. Grampp/Jaenicke
  1991 is a family-wide theory-vs-experiment compilation.
- **Couple:** Q/Q•⁻ (neutral/anion) — matches our reduction couples (anthraquinone 0/−1 etc.).
- **How to compare:** experimental total λ ≈ our λ_i(Q/Q•⁻) **+ λ_o(1e)**. For duroquinone,
  D3TaLES λ_i(electron)=0.573 eV + λ_o≈0.77 → predicted total ≈ 1.34 eV; check against
  Williams' k_ex-derived λ. Duroquinone is the ideal round-trip test because the SAME molecule
  has both a computed λ_i (D3TaLES) and an experimental self-exchange rate.
- **Caveats:** DDQ is strongly EWG-substituted — its λ ≠ parent BQ; use it as a
  substituent-effect anchor, not a parent proxy. Avoid PCET/protic couples (Q/QH₂): our λ_i is
  aprotic single-ET only. Solvent in Williams 1969 must be checked (may not be MeCN).

### TEMPO / nitroxides — ONE gold experimental source, and an interesting discrepancy to test
- **Direct self-exchange:** Grampp/Rasmussen 2002 (PCCP) — **TEMPO•/TEMPO+ ESR
  line-broadening self-exchange in a solvent series including MeCN**, with solvent-dynamical
  Marcus analysis. This is the oxidation (hole) couple.
- **How to compare:** our λ_i(TEMPO ox→rad)=0.971 eV and D3TaLES hole=0.943 eV **agree tightly
  (~0.03 eV)** — good internal consistency. Add λ_o≈0.78 → predicted total ≈ 1.75 eV.
- **FLAG / open question:** that computed inner-sphere λ_i (~0.95 eV) is **large**, yet TEMPO
  is famous for *fast, reversible* electrode kinetics. Either much of the geometric change is
  low-frequency/entropic (weakly rate-limiting) or the effective λ that sets k_ex is smaller
  than the vertical 4-point λ_i. **This is the single most valuable experimental check in the
  set** — pull Grampp/Rasmussen's derived total λ and reconcile. The TEMPO **electron** couple
  (rad/anion, λ_i=1.48 eV in D3TaLES) is almost certainly an **unbound/poorly-bound anion**
  artifact (TEMPO does not reduce cleanly); do not treat it as a real anchor.
- Mitra/Heuer/Diddens 2024 (PCCP) gives modern DFT+MD computed λ across solvents — a
  methodological cross-check for the solvent dependence, not an experimental anchor.

### Phenothiazines — gold IV-CT source, oxidation couple only
- **IV-CT (optical):** Sun/Rosokha/Kochi 2004 (JACS) — the intermolecular self-exchange
  precursor complex (PH)₂•⁺ and intramolecular bridged mixed-valence P(br)P•⁺; λ from the
  intervalence (charge-resonance) band by Mulliken–Hush two-state analysis. Rosokha/Kochi 2008
  (Acc. Chem. Res.) tabulates consolidated λ / H_ab for many organic donor couples — the best
  single place to pull numbers once you have access.
- **Couple:** PH/PH•⁺ (neutral/cation) — matches our phenothiazine **hole** couple.
- **How to compare:** our λ_i(phenothiazine ox→neu)=0.518 eV ≈ D3TaLES N-ethyl-PTZ hole
  0.45–0.63 eV (rep. 0.536 eV) — consistent. **Caveat:** IV-CT λ from a *contact* mixed-valence
  ion is dominated by inner-sphere + a *reduced* outer-sphere (short donor–acceptor distance),
  so an IV-CT λ is closer to our λ_i than a diffusive self-exchange total λ is. Match distance
  regime carefully; solvent in the IV-CT work is typically CH₂Cl₂, not MeCN.

### Viologens — WEAKEST experimental footing (flagged)
- No clean single-molecule primary DOI for MV self-exchange λ was pinned this session.
  Qualitative literature consensus (docs/LITERATURE_RFB.md, refs V1–V7): **MV2+/MV+• has low
  inner-sphere λ** (delocalized, planar radical cation → fast reversible ET); **MV+•/MV0 is
  higher** (planarization). Our pipeline: λ_i(2+/+•)=0.540, λ_i(+•/0)=0.373 eV. No D3TaLES
  cross-check (charged bipyridinium absent from the public dump).
- **Action needed:** targeted retrieval of (a) pulse-radiolysis / ESR self-exchange rate of
  MV•⁺, and (b) cyclophane/bridged **bis-viologen mixed-valence IV-CT** studies (the natural
  IV-CT anchor for this family). Until then, treat viologen λ as computed-only.
- **Calibration surrogate:** TMPD (Wurster's blue), Grampp/Jaenicke 1985, is the canonical
  low-λ aromatic-amine self-exchange and validates the workflow end-to-end even though it is
  not a target family.

---

## 3. Solid vs shaky — one-line summary

| Anchor | Status |
|---|---|
| D3TaLES computed λ_i (AQ, NQ, duroquinone, N-Et-PTZ, TEMPO hole) | **Solid** (real values, level of theory stated, same 4-point definition as ours) |
| Our pipeline λ_i (all families) | **Solid** as computed values; validated internally vs D3TaLES where SMILES overlap |
| Duroquinone: computed λ_i + Williams expt self-exchange | **Solid pairing** — best single validation target once k_ex is pulled |
| TEMPO Grampp/Rasmussen self-exchange (MeCN) | **Solid source, number to retrieve**; flagged discrepancy vs large computed λ_i |
| Phenothiazine Rosokha/Kochi IV-CT | **Solid source, number to retrieve**; distance/solvent regime caveat |
| DDQ, benzoquinone self-exchange | **Solid sources, numbers to retrieve**; DDQ ≠ parent |
| Viologen experimental λ | **Shaky/unretrieved** — no confident DOI; computed-only for now |
| Any numeric experimental λ in this CSV | **None filled** — all paywalled this session; retrieve from the cited DOIs |

## 4. Recommended level of theory to match the computed anchors

- To reproduce/extend the **D3TaLES** anchors exactly: **IP-tuned LC-ωHPBE / def2-SVP,
  gas-phase, 4-point Nelsen**, per-molecule tuned ω (Gaussian16). This is the level the
  in-repo D3TaLES cross-check already uses (`validate_reorg_worker_d3tales.py`); it removes the
  ~0.24 eV functional+basis offset seen vs our wB97M-V/def2-TZVP.
- Our production λ_i (**wB97M-V/def2-TZVP** on SMD geometries) is internally consistent and
  agrees with D3TaLES to ≈0.03–0.08 eV for AQ and TEMPO; keep it as the primary, and quote
  D3TaLES-level λ only when directly comparing to the D3TaLES database.
- To compare to **experimental total λ** (self-exchange / heterogeneous rate fits): report
  **λ_total = λ_i + λ_o** with λ_o from `solvated_reorg.py` in MeCN, and state explicitly
  whether the experiment was homogeneous self-exchange (use full λ_o) or heterogeneous
  (halve λ_o). Match the experiment's solvent (several classic studies are in CH₂Cl₂/DMF, not
  MeCN — the Pekar factor differs, so λ_o must be recomputed for that solvent).
- For **IV-CT** anchors (phenothiazine, bis-viologen): compare to λ_i more directly than to a
  diffusive total λ, because the contact mixed-valence ion has a compressed outer-sphere term.

## 5. Next actions to harden this set
1. Retrieve the numeric λ / k_ex from the 5 verified paywalled DOIs (institutional access):
   `10.1039/b206313a` (TEMPO), `10.1039/a903394g` (DDQ), `10.1080/00268976900100071`
   (BQ/duroquinone), `10.1021/ja038746v` + `10.1021/ar700256a` (phenothiazine).
2. Find and add a viologen self-exchange / bis-viologen IV-CT primary source.
3. Once numbers are in: for each, compute our λ_total in the *matching* solvent and tabulate
   (expt total) vs (λ_i + λ_o), reporting the residual as the validation metric.
