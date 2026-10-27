# Proposed phenothiazine linker and substitution matrix

[Structures](../../figures/candidates/phenothiazine_linker_matrix.png) · [Vector figure](../../figures/candidates/phenothiazine_linker_matrix.svg) · [Mapped designs](matrix.json) · [Route screen](route_screen.json) · [All-axis decisions](all_axis_decisions.json) · [Summary](summary.json) · [Tests](tests_verified.xml) · [Ledger check](ledger_identity_check.json)

This 2×2 hypothesis panel holds the phenothiazine O-linked local attachment class fixed while comparing a flexible two-carbon linker with one explicitly stereospecified cyclic linker, and unsubstituted rings with 3,7-dimethoxy substitution. The flexible unsubstituted product reuses the earlier [source-identified precursor proposal](../ptz_hydroxyethyl_source_2026-10-26/README.md); the other three graphs are new relative to that proposal and the frozen 879-graph ledger. Removing each linker gives the same matched phenothiazine core for a given substitution, and the two cyclic cases use the same assigned stereochemical pattern.

Methoxy substitution changes electronics as well as sterics; it is not evidence of shielding or reduced coupling. Ring closure changes topology and mass, and equivalent molecular reach or polymer conformer behavior has not been established. The flexible H/dimethoxy products are 347.48/407.54 g mol⁻¹ and the cyclic H/dimethoxy products are 359.49/419.55 g mol⁻¹. The cyclic cases specify one stereoisomer; precursor synthesis and stereochemical access are unresolved.

All four routes are **hold**, and all four shared all-axis material decisions are **pending**. The source-reported hydroxyethyl phenothiazine supports only the unsubstituted flexible precursor; the [3,7-dimethoxy phenothiazine study](https://doi.org/10.1039/C7TA05883G) concerns related parent molecules, not these exact grafts. Actual O-grafting conditions, redox retention, potential, grafting free energy, motif association, coupling and accessibility are unverified. One oxidation is proposed as a computational hypothesis; no reversibility, capacity or second electron is assumed. No quantum jobs were submitted, and no synthesis or material shortlist follows.

This is a single complementary oxidation-active comparison, not completion of all predeclared design slots or the broader benchmark. The three new products remain outside the previously frozen decision ledger.

Reproduce with fresh output names from the project root:

```bash
export PYTHONPATH="$PWD/src"
python scripts/pipeline/prepare_phenothiazine_linker_matrix.py --spec config/phenothiazine_linker_matrix_2026-10-26.json --baseline results/discovery/ptz_hydroxyethyl_source_2026-10-26/candidates.json --ledger results/selection/campaigns/20261024_UMA_default/decisions.json --out NEW_MATRIX.json --candidates-out NEW_CANDIDATES.json
python scripts/pipeline/screen_grafting_routes.py --candidates NEW_CANDIDATES.json --out NEW_SCREEN.json
python scripts/pipeline/assess_matched_route_panel.py --matrix NEW_MATRIX.json --gate NEW_SCREEN.json --policy config/candidate_feasibility.json --out NEW_DECISIONS.json
```

The [plotting script](../../../scripts/plotting/candidates/plot_phenothiazine_linker_matrix.py) draws the stored exact product graphs.

The route registry was later extended with the [anthraquinone carboxylate matrix](../anthraquinone_carboxylate_matrix_2026-10-26/README.md). Its [append-only reconciliation](../anthraquinone_carboxylate_matrix_2026-10-26/registry_reconciliation.json) preserves this matrix and the earlier registry entries.
