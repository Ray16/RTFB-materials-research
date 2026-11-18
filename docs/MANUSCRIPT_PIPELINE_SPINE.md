# Manuscript spine — pipeline, stage by stage (draft 2026-11-17)

Non-aqueous redox-flow-battery redox motifs grafted onto Merrifield (chloromethyl polystyrene)
polymer. Each pipeline stage is mapped to its manuscript role: the CLAIM we can make, the
EVIDENCE/figure that backs it, and the honest CAVEAT. Numbers are from the 92-candidate pool on
the uniform legacy xTB-RRHO thermal basis (MeCN, PF6-, vs Fc/Fc+).

## One-line thesis (to be confirmed)
A validated DFT+SMD screening pipeline grafts anolyte/catholyte motifs onto a Merrifield proxy
and ranks them on device axes; applied to 92 candidates it fills the missing 2e-catholyte side
with grafted dihydrophenazines carrying DFT-confirmed substituent design rules, and proposes a
~1.8 V 2e/2e cell.

---

## Stage 1 — Motif registration + grafting  (Methods; Fig. 1 schematic)
- CLAIM: motifs are grafted onto a common 4-methylbenzyl tether (proxy for the Merrifield
  benzylic displacement), redox states defined with explicit charge/multiplicity.
- EVIDENCE: config/custom_candidates.py (92 entries); molzip graft; SCAFFOLD Cc1ccc(C[*:1])cc1.
- CAVEAT: the tether is a proxy; exact polymer-bound synthesis/stability is out of scope (route
  review noted per candidate).

## Stage 2 — Conformer gen + UMA pre-optimization  (Methods)
- CLAIM: MLIP (UMA) pre-opt over conformers controls cost before DFT.
- CAVEAT: cosmetics only; no property computed here.

## Stage 3-4 — DFT geometry + SMD, single-points  (Methods: level of theory)
- CLAIM: wB97M-V/def2-TZVPD with SMD(MeCN), density-fitted, per redox state.
- EVIDENCE: archive calcs/dft/<id>/<state>/result.json (e_smd_eV, e_gas_eV).

## Stage 5 — Thermochemistry  (Methods caveat)
- CLAIM: Gibbs energies on a uniform legacy xTB-RRHO basis across all 92.
- EVIDENCE: backends.dft._thermal_correction (gfn2 RRHO).
- CAVEAT: native solution-DFT qRRHO is implemented but NOT run pool-wide (compute/disk). Shown:
  candidate-side native-vs-legacy thermal is identical (0 mV) -> comparisons valid on one basis.

## Stage 6 — Properties  (CORE RESULTS)
| property | coverage | reliability | figure/table |
|----------|----------|-------------|--------------|
| redox potential (vs Fc) | 92/92 | MAE 0.226 V; non-uniform bias | validation parity |
| usable capacity (TBAPF6/MeCN) | 92/92 | exact arithmetic | landscape fig |
| outer-sphere lambda_o (Born) | 92/92 | continuum, shape-blind | lambda_o table |
| inner-sphere lambda_i (4-point) | catholytes + partial | reliable for OXIDATIONS; low-conf for reductions (anion_unbound) | reorg table |
| disproportionation stability | per couple | from E1,E2 spacing | — |

## Stage 7 — QC gating  (Methods: reliability)
- CLAIM: spin, geometry integrity, state coverage, anion-unbound flags gate what is scored.
- CAVEAT: anion_unbound is a real physical flag (n-type reductions), not a bug; it is why
  lambda_i is held for anolytes.

## Stage 8 — Validation  (CREDIBILITY ANCHOR; parity figure)
- CLAIM: potentials reliable to MAE 0.226 V with a correctable, STRATIFIED bias.
  - overall bias +0.191 V; reductions +0.166, oxidations +0.268; by family imide +0.072 (best),
    quinone +0.253 (worst).
- CATHOLYTE ANCHOR: computed DHP oxidations sit +0.30/+0.44 V above measured BMEPZ (-0.29/+0.50);
  consistent with the oxidation bias + a real substituent shift -> catholyte potentials biased HIGH.
- CAVEAT (reviewer-facing): no same-chemistry (phenazine/amine) couple in the tier-A benchmark, so
  the catholyte bias is CORROBORATED, not a per-couple MAE. The clean experiment to close this is
  computing the exact BMEPZ parent (queued; disk-gated).

## Stage 9 — Scoring + Pareto proposal  (PROPOSAL FIGURE)
- CLAIM: sigma-aware Pareto fronts on RAW axes (potential x usable capacity), no scalarized FoM.
  Anolyte front 7, catholyte front 3; the expansion fills the catholyte side.
- KEY RESULT (the expansion's contribution): grafted dihydrophenazine 2e catholytes with
  DFT-confirmed design rules:
  - smaller N'-cap -> higher capacity (N'-methyl-DHP 172 mAh/g active repeat-unit mass, pool best 2e catholyte)
  - EWG N'-cap -> higher potential (N'-CF3CH2-DHP +0.62 V)
  - ring donor -> lower potential, MONOTONIC: unsub +0.43 -> 2,7-diMe +0.30 -> 2,7-diOMe +0.21 V
- PROPOSED LEAD CELL: 2-amino-1,4-NQ (anolyte -1.41 V, 185 mAh/g active, lambda_o 0.73 eV) +
  N'-methyl-DHP (catholyte +0.47 V, 172 mAh/g active, lambda_o 0.70 eV, lambda_i reliable)
  -> ~1.8 V (1.78-1.87 bias-corrected), +/- ~0.4 V per-couple uncertainty.
- CAVEAT: a full CELL has no experimental datum (candidates are un-synthesized proposals); only
  the component potentials are benchmarked.

---

## Figures/data assets already produced
- results/figures/candidates/candidates_potential_vs_capacity_MeCN_TBAPF6.png  (landscape + fronts)
- results/validation/full_pool_scored_legacyG_MeCN_TBAPF6.csv  (92 rows)
- results/validation/lambda_outer_fullpool_2026-11-17.csv  (92 rows)
- results/validation/catholyte_reorg_and_anchor_2026-11-17.md  (reorg + bias + anchor)
- archive/0_redox/results/validation/benchmark_mecn.csv  (44 tier-A couples)

## Open items (all gated on freeing the 100%-full shared NFS, or impossible by construction)
1. Native qRRHO thermochemistry pool-wide  (compute+disk).
2. Exact BMEPZ parent -> per-couple catholyte MAE  (compute+disk).
3. Experimental full-cell validation  (impossible: candidates un-synthesized; validate components).
