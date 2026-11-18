> **Superseded (2026-11-18).** This per-section draft is kept for history only and is out of date:
> several numbers and claims in it were later corrected. The current text is
> `manuscript/manuscript_draft_v0.md` (journal version: `manuscript/jcim/`; SI:
> `manuscript/supporting_information.md`). Do not cite or copy from this file.

## Screening the candidate space, the catholyte gap, and proposed full cells

### Mapping 91 grafted candidates onto a two-pool decision front

Applied to the full set of 91 registered grafted motifs, the pipeline returns a
reliability-gated property record for every candidate — potential versus
ferrocene (Fc) in TBAPF₆/MeCN, nominal usable capacity, outer- and inner-sphere
reorganization estimates, and disproportionation free energies — from which the
accessible design space can be read directly (Figure 5). The pool is strongly
anolyte-dominated: 85 of the 91 candidates reduce (n-type anolytes) and only 6
oxidize into an accessible catholyte window (p-type). By family the set is 36
imides, 29 quinones, 18 multi-electron pyridines, 5 dihydrophenazines, 2
mono-imides, and 1 phenothiazine; 89 of the 91 are two-electron couples, the two
exceptions being the single-wave phthalimide and 1,8-naphthalimide controls.
This composition is itself a finding about the literature we draw from:
graftable, well-precedented n-type cores (naphthalene- and pyromellitic-diimides,
benzo-/naphtho-/anthraquinones) are abundant, whereas graftable p-type cores are
scarce.

Candidates are ranked not by a single figure of merit but on σ-aware Pareto
fronts computed **separately for the anolyte and catholyte pools**, because a
redox potential is not globally "more is better": an anolyte is better the more
negative its potential, a catholyte the more positive. Within each pool the
trustworthy axes — pool-appropriate voltage and nominal usable capacity (mAh g⁻¹
in TBAPF₆/MeCN) — define the objectives, oriented so that higher is better, and
domination is declared only when one candidate is no worse than another beyond
the combined per-candidate uncertainty (σ, propagated from the validation
benchmark) on every shared axis and strictly better, beyond that noise, on at
least one. Candidates missing a primary axis are reported on a separate
partial-information list rather than silently promoted or dropped. We
deliberately avoid any weighted-sum score: arbitrary weights, insensitivity to
σ, and pool-relative min–max rescaling would all make a scalar rank change
whenever the candidate set changes. The reader is given the front and the raw
axes instead (Figure 5; full table in the SI). Inner-sphere reorganization
(λᵢ) is carried only qualitatively and is not used to set the headline ranking,
for the reasons given in the Methods.

Figure 5 presents the two Pareto panels side by side. The anolyte panel (85
candidates) is densely populated across a wide potential window: diimides and
amino-/alkoxy-naphthoquinones span roughly −1.0 to −1.6 V vs Fc, with
pyromellitic diimides (e.g. the ether-tethered PMDI at −1.54 V) and
alkoxy-anthraquinones (≈ −1.49 V) anchoring the most-negative edge and the
multi-electron viologens clustering nearer −0.7 to −1.0 V. Nominal capacities
across the anolyte pool fall in the ≈ 50–72 mAh g⁻¹ band, trending down as
solubilizing oligoether or long-alkyl caps add inert mass. The catholyte panel,
by contrast, is almost empty: six points, five of them dihydrophenazines
clustered near +0.4 to +0.5 V and one 3,7-dimethoxyphenothiazine at +0.78 V.
This visual sparsity — a crowded anolyte front against a nearly bare catholyte
front — is the central structural result of the screen, and it defines the
**catholyte gap** that the next step is designed to fill.

### Filling the catholyte gap with grafted dihydrophenazines

To populate the empty catholyte front we used the pipeline's matched-expansion
machinery to graft a single scaffold — 5,10-dihydrophenazine (DMPZ) — and vary
one substituent at a time, so that the resulting potential shifts can be read as
matched single-substituent deltas rather than confounded whole-molecule
comparisons (Figure 6). Every member carries the same 4-methylbenzyl (Merrifield
proxy) cap on one nitrogen and a two-electron {neutral → +1 → +2} ladder,
isolating the effect of the N′-substituent and of ring substitution. The
confirmed series potentials (E vs Fc, mean of the two one-electron couples, from
`full_pool_scored_legacyG_MeCN_TBAPF6.csv`) are:

| Candidate | N′ / ring substitution | E vs Fc (V) | Couples (V) | Capacity (mAh g⁻¹) |
|---|---|---:|---|---:|
| `exp_dmpz_mb_ethyl` | N′-ethyl | +0.43 | −0.03 / +0.89 | 88.7 |
| `exp_dmpz_mb_methyl` | N′-methyl | +0.47 | +0.03 / +0.91 | 90.8 |
| `hyp_dmpz_mb_methoxyethyl` | N′-(2-methoxyethyl) | +0.47 | +0.01 / +0.94 | 84.5 |
| `hyp_dmpz_mb_trifluoroethyl` | N′-(2,2,2-trifluoroethyl) | +0.62 | +0.15 / +1.10 | 81.4 |
| `hyp_dmpz_27dimethyl_ethyl` | 2,7-dimethyl, N′-ethyl | +0.30 | −0.14 / +0.73 | 84.8 |
| `hyp_dmpz_27dimethoxy_ethyl` | 2,7-dimethoxy, N′-ethyl | +0.21 | −0.22 / +0.64 | 80.7 |
| `exp_ptz_37dimethoxy` | 3,7-dimethoxyphenothiazine | +0.78 | +0.37 / +1.19 | 82.0 |

Read as matched deltas against the N′-ethyl DMPZ reference (+0.43 V), the series
yields a clean, transferable **design rule**. An electron-withdrawing
N′-substituent raises the potential: replacing ethyl with 2,2,2-trifluoroethyl
shifts E by +0.20 V (to +0.62 V), while the near-isosteric N′-methyl and
N′-(2-methoxyethyl) caps barely move it (+0.04 to +0.05 V). Conversely,
electron-donating ring substitution lowers it: 2,7-dimethyl brings E down by
≈ 0.13 V (to +0.30 V) and 2,7-dimethoxy by ≈ 0.22 V (to +0.21 V) relative to the
same N′-ethyl parent. The two handles are roughly additive and act on opposite
ends of the scaffold, so N′-electron-withdrawal and ring-electron-donation give
independent, oppositely-signed control over where a DMPZ catholyte sits — spanning
≈ +0.21 to +0.62 V across the grafted series, with the distinct
3,7-dimethoxyphenothiazine scaffold extending the catholyte reach to +0.78 V.
Equally important for a flow-battery catholyte, these light two-electron
phenazine cores deliver high nominal capacity (≈ 81–91 mAh g⁻¹), above the
anolyte band, because the redox-active mass fraction is large. Re-plotting the
catholyte Pareto panel with these points shows the gap being filled by a
structured family rather than by isolated outliers (Figure 6).

Several of these values are provisional. The expansion campaign
(`results/expand_catholyte_20261117`) is live at the time of writing, with the
N′-methyl, N′-(2,2,2-trifluoroethyl), 2,7-dimethyl, and 2,7-dimethoxy DMPZ
members still completing their native-DFT confirmation stages; their legacy
ranking potentials are reported here and will be reconciled against the
disk-gated qRRHO values once the run finishes and QC gates pass. The frozen
91-candidate pool underlying Figure 5 contains six catholytes (five DMPZ + one
phenothiazine); the maintained scored CSV now carries a seventh catholyte — the
2,7-dimethoxy DMPZ (E ≈ +0.21 V) appended by the live campaign — so it holds 92
rows while the frozen analysis pool remains 91. The figures enforce the frozen 91
and footnote the live extra.

### Proposed non-aqueous full cells

Because every potential in both pools is referenced to the same validated Fc
scale, candidates can be paired directly and a nominal cell voltage read off as
E(catholyte) − E(anolyte) (Figure 7). Pairing a surfaced DMPZ catholyte
(≈ +0.4 to +0.6 V) with a pipeline anolyte from the most-negative edge of the
quinone/imide front (≈ −1.5 V) gives proposed non-aqueous full cells near
**~1.9 V** — for example N′-ethyl or N′-methyl DMPZ (+0.43 / +0.47 V) against an
alkoxy-anthraquinone or pyromellitic-diimide anolyte (≈ −1.48 V) yields ≈ 1.9–1.95 V,
and the more electron-poor N′-(2,2,2-trifluoroethyl) DMPZ (E ≈ +0.62 V) paired
with an ether-tethered PMDI anolyte (≈ −1.48 to −1.54 V) would extend the proposal
toward ≈ 2.1 V (mean-potential arithmetic verified against the scored CSV:
+0.43/+0.47 − (−1.48) = 1.91/1.95 V; +0.62 − (−1.48) = 2.10 V). The design rule above makes
this voltage a tunable quantity rather than a single number: catholyte and
anolyte potentials can be shifted independently to target a window.

We state the limitation plainly. A computed cell voltage derived from two
validated half-cell potentials is a **proposal**, not a demonstrated device. It
says nothing about attainable capacity under cycling, Coulombic efficiency,
membrane crossover of the small grafted proxies, chemical stability of the
charged states over many cycles, or the synthetic feasibility of the exact
grafted structure — each candidate's record carries explicit route caveats, and
a capped graft is not its soluble parent nor an oligomer. What the screen does
establish is a reliability-gated, internally consistent map of where these
chemistries sit, a transferable substituent rule for moving a catholyte to a
target potential, and a concrete, testable set of anolyte/catholyte pairings for
experimental follow-up.
