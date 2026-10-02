# two_wave_mecn — experimental 1st/2nd one-electron potentials in MeCN

`two_wave_mecn.csv`: literature E1 (first one-electron wave) and E2 (second wave) for quinones,
aromatic imides and viologens, measured in **acetonitrile with a tetraalkylammonium electrolyte**.
Curated 2026-10 from primary papers that were actually opened (open-access full text, author
manuscripts, or ACS SI files mirrored on acs.figshare.com). 37 rows: 21 have both waves, 3 are
one-wave phthalimides, and 13 are first-wave-only quinones from one internally consistent set
(Huynh et al. 2016).

## Columns
`E*_orig_V` and `ref_orig` hold the value and reference **exactly as the paper reports them**.
`E*_vs_Fc_V` holds the value on the Fc+/Fc scale in MeCN. `dE12_V = E1 − E2` does not depend on the
reference. For viologens, `smiles_neutral` holds the **dication** (this is said in `notes`).
`confidence`: high = E1/2 or E°' measured directly against internal Fc; medium = peak or half-peak
potentials, or a conversion through a reference electrode; low = identity, solvent or wave is ambiguous.

## Conversion rules (applied in this order)
1. If the paper reports values vs Fc+/Fc (internal ferrocene), they are used as-is.
2. If the paper converted from its own measured Fc onto another scale, the authors' own constant
   is inverted exactly:
   - Schimanofsky 2022 used Fc = +0.620 V vs SHE in MeCN, so E_Fc = E_SHE − 0.620.
   - Cabral 2023 used Fc = +0.40 V vs SCE, so E_Fc = E_SCE − 0.40.
3. If no Fc value is given (only Cook & Horrocks 2017, methyl viologen, Ag/Ag+ 10 mM), the
   conversion is E_Fc = E_Ag/Ag+(0.01 M AgNO3) − 0.082 V. This comes from Pavlishchuk & Addison,
   Inorg. Chim. Acta 2000, 298, 97 (plus corrigendum doi 10.1016/j.ica.2024.122468). **SECONDARY:** the
   primary paper is paywalled and was not seen. The constant is taken from the Lam group converter
   (lamresearchgroup.com/potential-converter), which uses offsets Fc −37 mV, Ag/Ag+(0.01 M) +45 mV and
   SCE +343 mV relative to Ag/Ag+(0.1 M), for 0.1 M TEAP/Et4NPF6.
   - The same source gives **Fc = +0.380 V vs SCE**, not the +0.40 often attributed to P&A.
   - Primary measurements we did see: Fc = +0.399 V vs SCE (MeCN/0.1 M Et4NBF4; Aguilar-Martínez,
     Macías-Ruvalcaba, González, Rev. Soc. Quím. Méx. 2000, 44, 17). Cabral 2023 also used 0.40 V.
   - The often-quoted "−0.087 V" Ag/Ag+ → Fc constant was **not** verified.
   - Commercial Ag/Ag+ electrodes vary: Frontana 2005 measured Fc at +0.25 V vs its own
     Ag/0.01 M AgNO3 (TBAP) reference. Any Ag/Ag+ → Fc conversion therefore carries ~0.1 V risk.

## Sources used
| Data | Source | Where |
|---|---|---|
| BQ, 2-Me-BQ, 2-tBu-BQ, 2-F-BQ | Wang et al., J. Phys. Chem. C 2020, 124, 13609 | SI Table S1 |
| Cross-check (BQ) | Hooe et al., Chem. Sci. 2021, 12, 9733 | — |
| E1 only, 13 quinones, plus E1 cross-checks | Huynh et al., JACS 2016, 138, 15903 | SI Table S3 |
| 6 anthraquinones (Ep/2) | Schimanofsky et al., J. Phys. Chem. C 2022, 126, 14138 | Table 1 |
| 1,4-NQ and 1,2-NQ | Cabral et al., Molecules 2023, 28, 1232 | Table 3 |
| Juglone, naphthazarin, NQ cross-check (peak potentials) | Frontana & González, J. Mex. Chem. Soc. 2005, 49, 61 | Table 1 |
| PMDIs, phthalimides | Daub, Janssen, Hendriks, ACS Appl. Energy Mater. 2021, 4, 9248 | SI Table S1 |
| Methyl viologen | Cook & Horrocks, ChemElectroChem 2017, 4, 320 | Table 1 |
| 1-Benzyl-1'-methyl viologen | Kim, Sanford, Vaid, McNeil, Chem. Eur. J. 2022, 28, e202200149 | — |
| NDI (low confidence) | Xiao et al., Chem. Sci. 2022, 13, 13426 | — |

## Excluded (not MeCN/tetraalkylammonium, secondary, or ambiguous)
- **Prince, Dutton & Gunner, BBA Bioenerg. 2022 (350 quinones, two waves):** DMF vs SCE. Wrong
  solvent. It is the best large two-wave set if a DMF track is ever wanted.
- **Bui/Forse, J. Phys. Chem. C 2022 (AQ derivatives):** DMSO.
- **"Re-assessing viologens", Chem. Sci. 2024; Materials 2021, 14, 1702 (alkyl/benzyl viologens):**
  aqueous.
- **Chem. Sci. 2016 (PMC6013826), MV −0.89 V vs Fc:** secondary citation in a DMF study.
- **Angew. Chem. 2026 (PMC12851006), "paraquat −1.09/−1.52 V vs Fc":** cites Bird & Kuhn 1981
  and is consistent with an aqueous value pushed onto the Fc scale. Not used.
- **Martinez/La Porte/Wasielewski, J. Phys. Chem. C 2018 (NDI):** DMF.
- **Other excluded data:**
  - 2-hydroxy-1,4-NQ (self-protonation).
  - Wang 2020 di-substituted quinones (2,5 vs 2,6 isomer not determinable from the text).
  - Daub mellitic triimides (three waves).
  - ANDI (Li/K TFSI electrolyte).
- **Bird & Kuhn, Chem. Soc. Rev. 1981 (viologens):** paywalled, not seen. Its MV2+/MV+• value of
  −0.446 V vs NHE is the aqueous couple, so it must not be used as an MeCN anchor.

## Caveats
- The 1,4-naphthoquinone E2 disagrees between sources by 0.29 V (−1.47 vs −1.76 V vs Fc;
  Cabral's value is a cathodic peak). Treat it as unreliable.
- Schimanofsky values are Ep/2, which sit about 0.03 V positive of E1/2.
- Stereocenters: the 2-ethylhexyl derivatives contain unspecified stereocenters. Enumerate the
  stereoisomers before any calculation (project rule).

## Update (follow-up): oxidation anchors + new columns
- New columns: `couple` (charge states, e.g. "0/-1; -1/-2", "+2/+1; +1/0", "+1/0" for oxidations; for
  oxidations E1 = first oxidation potential) and `wave2_reversible` (true only if the source gives an
  E1/2/E°' for wave 2 from a reversible/quasi-reversible couple; false if wave 2 is only Epc/Ep/2,
  described as poorly behaved, or not reported at all - one-wave and E1-only rows are false). For
  Daub 2021 and Kim 2022, "true" is inferred from reported E1/2 values (reversibility not described
  in words).
- TEMPO+/TEMPO: Gerken & Stahl, ACS Cent. Sci. 2015 (+0.249 V vs Fc, reversible; electrolyte 0.5 M
  KPF6, NOT tetraalkylammonium) and Bingham et al., JACS 2026 (peaks 0.26/0.17 V vs Fc in
  0.1 M Bu4NClO4/MeCN with 0.8 M H2O -> midpoint 0.215 computed here). Both support the +0.24 anchor
  to within ~0.03 V. The config's citation "ACS Cent. Sci. 2015, 1, 234" is Gerken & Stahl.
- Phenothiazines: no clean primary MeCN/TAA E1/2 vs Fc found. Rows are low confidence:
  10-methylphenothiazine 0.75 V vs SCE (Mayther/Vullev JPCB 2023; an extrapolated neat-solvent value)
  -> 0.37 V vs Fc via secondary P&A constant; 0.34 V vs Fc for MePT and 0.23 V vs Fc for 10H-PTZ are
  secondary (Ishimatsu review 2025 quoting JEAC 2024, not seen). Both bracket the +0.26 anchor only
  loosely; the MePT anchor is therefore unverified.
- N-methylpyridinium (parent) reduction in MeCN: NO primary source found (open-access search found
  only substituted/extended pyridiniums). No row added; the −1.8 V anchor remains uncited.
- Pavlishchuk & Addison 2000 (and 2024 corrigendum): ScienceDirect/ResearchGate/ScienceOpen all
  returned 403; Crossref/Semantic Scholar have no abstract. Paper NOT accessed; constants remain
  secondary (Lam converter: Fc = SCE − 0.380 V; Fc = Ag/Ag+(0.01 M AgNO3) − 0.082 V).
- Additional conflict noted (not tabulated): Chem. Eur. J. 2021 (doi 10.1002/chem.202004748,
  PMC7898908) quotes paraquat −1.076/−1.51 V vs Fc "in MeCN/[Et4N][PF6]" via citations; this is
  ~0.23 V more negative than the direct MeCN measurements (Kim 2022, Cook 2017) and matches an
  aqueous-derived value, so it was not used.
