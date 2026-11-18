# Figure 6 — Grafted dihydrophenazine (DMPZ) catholyte design rule

**Figure 6. A matched single-substituent rule for tuning a grafted
5,10-dihydrophenazine (DMPZ) catholyte.** Mean redox potential versus ferrocene
(E vs Fc/Fc⁺) in TBAPF₆/MeCN for a family built by grafting one scaffold — DMPZ
with a 4-methylbenzyl (Merrifield-proxy) cap on one nitrogen and a two-electron
{neutral → +1 → +2} ladder — and varying a single substituent at a time, so that
each shift reads as a matched delta rather than a confounded whole-molecule
comparison. Potentials (E = mean of the two one-electron couples) are read live
from `results/validation/full_pool_scored_legacyG_MeCN_TBAPF6.csv`; no value is
altered. DMPZ members are drawn in orange and the distinct
3,7-dimethoxyphenothiazine (PTZ) in purple (fixed family colours across the
figure set). Members are ordered by potential; the dashed line marks the N′-ethyl
DMPZ reference (+0.43 V), and the parenthetical value in each annotation is the
matched shift relative to that reference (red where a substituent raises E, green
where it lowers E).

The series gives a simple design rule with two handles that push the potential in
opposite directions. An electron-**withdrawing**
N′-substituent **raises** E: N′-(2,2,2-trifluoroethyl) sits +0.20 V above the
N′-ethyl reference (to +0.62 V), while the near-isosteric N′-methyl and
N′-(2-methoxyethyl) caps barely move it (+0.04 to +0.05 V). An
electron-**donating** ring substituent **lowers** E: 2,7-dimethyl by −0.13 V (to
+0.30 V) and 2,7-dimethoxy by −0.22 V (to +0.21 V) relative to the same N′-ethyl
parent. The grafted DMPZ series therefore spans +0.21 to +0.62 V, and the
distinct PTZ scaffold reaches +0.78 V. Whether the two handles add when combined in
one molecule was not tested, since no member carries both. Panel **(b)**
re-plots the catholyte sub-space (E vs Fc/Fc⁺ against theoretical capacity on the
neutral repeat-unit mass); these light two-electron cores reach ~140–170 mAh g⁻¹,
comparable to the anolyte pool, and the structured family fills the catholyte gap rather than
scattering as isolated outliers.

All potentials use the pool's ranking thermal tier (GFN2-xTB RRHO), like every
value in Figure 5; none of the DMPZ members has yet been confirmed with the more
expensive native-Hessian tier. A check against the measured parent compound
(BMEPZ) suggests the absolute DMPZ potentials are about 0.3 V too high; the shifts
relative to the reference are the robust quantity. Potentials are QC-gated values
on the Fc/Fc⁺ scale (benchmark MAE 0.247 V over n = 39 endpoints; the error is
mostly a shared offset, so ranking is more reliable than absolute values); a
computed descriptor or a completed job establishes neither synthetic feasibility
nor material performance, and a capped graft is not its soluble parent nor an
oligomer.
