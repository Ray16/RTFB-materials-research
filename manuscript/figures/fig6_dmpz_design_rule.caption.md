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

The series yields a clean, roughly additive design rule acting through two
independent handles on opposite ends of the scaffold. An electron-**withdrawing**
N′-substituent **raises** E: N′-(2,2,2-trifluoroethyl) sits +0.20 V above the
N′-ethyl reference (to +0.62 V), while the near-isosteric N′-methyl and
N′-(2-methoxyethyl) caps barely move it (+0.04 to +0.05 V). An
electron-**donating** ring substituent **lowers** E: 2,7-dimethyl by −0.13 V (to
+0.30 V) and 2,7-dimethoxy by −0.22 V (to +0.21 V) relative to the same N′-ethyl
parent. The grafted DMPZ series therefore spans +0.21 to +0.62 V, and the
distinct PTZ scaffold extends the catholyte reach to +0.78 V. Panel **(b)**
re-plots the catholyte sub-space (E vs Fc/Fc⁺ against usable capacity in
TBAPF₆/MeCN); these light two-electron cores deliver high usable capacity
(~81–91 mAh g⁻¹), and the structured family fills the catholyte gap rather than
scattering as isolated outliers.

Open markers and an asterisk denote **provisional** members: the matched
expansion campaign (`results/expand_catholyte_20261117`) is live at the time of
writing and the native-DFT confirmation stages of the N′-methyl,
N′-(2-methoxyethyl), N′-(2,2,2-trifluoroethyl), 2,7-dimethyl, and 2,7-dimethoxy
DMPZ members are still completing. The legacy xTB-RRHO ranking potentials are
plotted here and will be reconciled against the disk-gated qRRHO values once the
run finishes and QC gates pass. Potentials are reliability-gated descriptors
referenced to the validated Fc/Fc⁺ scale (benchmark MAE 0.247 V over n = 39
endpoints; systematic bias so that ranking is more reliable than absolutes); a
computed descriptor or a completed job establishes neither synthetic feasibility
nor material performance, and a capped graft is not its soluble parent nor an
oligomer.
