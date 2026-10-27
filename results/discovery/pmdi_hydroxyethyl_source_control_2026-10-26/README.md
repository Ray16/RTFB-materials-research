# Source-reported hydroxyethyl PMDI precursor

[Structure and competing handles](../../figures/candidates/pmdi_hydroxyethyl_source_control.png) · [Vector figure](../../figures/candidates/pmdi_hydroxyethyl_source_control.svg) · [Mapped candidate](candidates.json) · [Source review](source_review.json) · [Route gate](route_screen.json) · [Material decision](all_axis_decision.json) · [Checks](summary.json) · [Tests](tests_verified.xml)

[US4320129A, Preparation 4](https://patents.google.com/patent/US4320129A/en) reports isolation of N-(2-hydroxyethyl)pyromellitic diimide. We interpret that name as the depicted mono(hydroxyethyl) PMDI with one alcohol and one imide N–H. The graph has formula C₁₂H₈N₂O₅ and is distinct from the N-ethyl-capped precursor in the [linker matrix](../pmdi_linker_matrix_2026-10-26/README.md). The source is one patent report; it does not establish current procurement or independent replication.

The proposed capped-site product attaches the alcohol oxygen to a benzyl carbon and retains the imide N–H. The mapped neutral-reference reaction is atom and charge balanced, but the alcohol and N–H provide competing reactive handles. The related [TEMPO bead attachment](../../grafting/tempo_exact_support_2026-10-26/README.md) concerns a different precursor and does not demonstrate selectivity or conditions for this PMDI compound.

The practical route screen remains **hold** and the all-axis decision is **pending**. Exact O-selectivity, conditions, redox retention, reaction thermodynamics, properties, interaction and accessibility are unresolved. No quantum jobs were submitted. This exact product is absent from the previously frozen 879-graph ledger and has not been inserted into it.

Reproduce with fresh output names: `PYTHONPATH=src python scripts/pipeline/prepare_pmdi_hydroxyethyl_control.py --out NEW_CANDIDATE.json`, followed by `scripts/pipeline/screen_grafting_routes.py` and `scripts/pipeline/assess_pmdi_hydroxyethyl_control.py`. The [plotting script](../../../scripts/plotting/candidates/plot_pmdi_hydroxyethyl_control.py) draws the mapped precursor and product.
