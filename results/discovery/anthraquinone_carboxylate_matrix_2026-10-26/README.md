# Proposed anthraquinone carboxylate matrix

[Structures](../../figures/candidates/anthraquinone_carboxylate_matrix.png) · [Vector figure](../../figures/candidates/anthraquinone_carboxylate_matrix.svg) · [Mapped designs](matrix.json) · [Route screen](route_screen.json) · [All-axis decisions](all_axis_decisions.json) · [Summary](summary.json) · [Tests](tests_verified.xml) · [Ledger check](ledger_identity_check.json)

This 2×2 hypothesis panel keeps the proposed benzyl-ester attachment class while comparing direct carboxylate attachment with one explicitly stereospecified cyclic spacer, and an unsubstituted anthraquinone ring with an explicit methyl substitution. The direct unsubstituted graph reuses the earlier [literature-linked AQ-2-carboxylate proposal](../aq_carboxylate_2026-10-26/README.md); the other three are new relative to that proposal and the frozen 879-graph ledger. Removing the attachment portion gives the same anthraquinone core for a given substitution, and the two cyclic cases use the same assigned stereochemical pattern. Both quinone carbonyls are retained in every graph.

The spacer changes topology and interrupts the direct carboxyl-to-core conjugation; methyl also changes electronics and sterics. Thus this is a design comparison, not an isolated measurement of linker rigidity or shielding. Equivalent reach, precursor synthesis, stereochemical access, ester stability and polymer conformer behavior are unverified. The direct H/Me products are 356.38/370.40 g mol⁻¹ and the cyclic H/Me products are 396.44/410.47 g mol⁻¹.

All four routes are **hold**, and all four shared all-axis material decisions are **pending**. [Oyaizu and coauthors](https://doi.org/10.1002/pat.1968) report attachment of the baseline acid to soluble chloromethylated polymer, not the three new exact precursors or crosslinked beads. Actual grafting conditions, selectivity, redox retention, potentials, grafting free energy, motif association, coupling and accessibility require further evidence. One reduction is proposed as a computational hypothesis; no reversibility or usable capacity is assumed. No quantum jobs were submitted and no synthesis or material shortlist follows.

This is a single AQ comparison, not completion of all predeclared slots or the broader benchmark. The three new products remain outside the frozen decision ledger.

Reproduce with fresh output names from the project root:

```bash
export PYTHONPATH="$PWD/src"
python scripts/pipeline/prepare_anthraquinone_carboxylate_matrix.py --spec config/anthraquinone_carboxylate_matrix_2026-10-26.json --baseline results/discovery/aq_carboxylate_2026-10-26/candidates_mapped.json --ledger results/selection/campaigns/20261024_UMA_default/decisions.json --out NEW_MATRIX.json --candidates-out NEW_CANDIDATES.json
python scripts/pipeline/screen_grafting_routes.py --candidates NEW_CANDIDATES.json --out NEW_SCREEN.json
python scripts/pipeline/assess_matched_route_panel.py --matrix NEW_MATRIX.json --gate NEW_SCREEN.json --policy config/candidate_feasibility.json --out NEW_DECISIONS.json
```

The [plotting script](../../../scripts/plotting/candidates/plot_anthraquinone_carboxylate_matrix.py) draws the stored exact product graphs.
