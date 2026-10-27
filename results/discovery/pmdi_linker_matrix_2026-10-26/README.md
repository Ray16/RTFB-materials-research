# Proposed PMDI linker and edge-substituent matrix

[Structures](../../figures/candidates/pmdi_linker_matrix.png) · [Vector figure](../../figures/candidates/pmdi_linker_matrix.svg) · [Mapped designs](matrix.json) · [Route screen](route_screen.json) · [All-axis decisions](all_axis_decisions.json) · [Summary](summary.json) · [Tests](tests_verified.xml) · [Ledger check](ledger_identity_check.json)

This 2×2 hypothesis panel fixes the ethyl-capped PMDI core and O-linked local attachment while varying a flexible two-carbon linker versus one explicitly stereospecified cyclopropyl linker, and H versus an edge methyl group. The two flexible graphs already occur in the systematic library; the cyclic graphs are two new exact proposals. The chosen cyclic stereoisomer has two assigned stereocenters, and its precursor synthesis and stereochemical access are not established.

The methyl group is a proposed steric perturbation, not evidence that stacking or electronic coupling is reduced. The two linker types change topology and mass: the flexible H/Me products are 392.41/406.44 g mol⁻¹ and the cyclic H/Me products are 404.42/418.45 g mol⁻¹. RDKit reports six versus five rotatable bonds, respectively; this descriptor does not establish equivalent reach, actual conformational ensembles or polymer behavior. Removing the linker from each graph gives the same matched PMDI core for a given edge substitution.

All four prospective route screens are **hold**, and the shared all-axis material policy returns **pending** for each. The cited patent reports a related monoethyl PMDI precursor, and the phthalimide polymer example supports a different reaction class; neither demonstrates these exact hydroxy precursors or their bead grafting. Exact precursor procurement, stereoselectivity, conditions, core retention, actual-route ΔG, redox properties, motif interactions and accessibility remain unresolved. The neutral HCl reactions are balanced reference bookkeeping, not established actual-condition reactions. No quantum jobs were submitted and no synthesis or material shortlist follows.

This is a concrete candidate panel toward the predeclared expansion, not completion of all 34 design slots or the broader benchmark. The new graphs are kept outside the previously frozen 879-graph decision ledger, which is not silently changed.

Reproduce with fresh output names from the project root:

```bash
export PYTHONPATH="$PWD/src"
python scripts/pipeline/prepare_pmdi_linker_matrix.py --library results/discovery/systematic_families_2026-10-18/structures.json --out NEW_MATRIX.json --candidates-out NEW_CANDIDATES.json
python scripts/pipeline/screen_grafting_routes.py --candidates NEW_CANDIDATES.json --out NEW_SCREEN.json
python scripts/pipeline/assess_pmdi_linker_matrix.py --matrix NEW_MATRIX.json --gate NEW_SCREEN.json --policy config/candidate_feasibility.json --out NEW_DECISIONS.json
```

The [plotting script](../../../scripts/plotting/candidates/plot_pmdi_linker_matrix.py) draws the stored exact product graphs.

The route registry was later extended with the [hydroxyethyl PMDI source control](../pmdi_hydroxyethyl_source_control_2026-10-26/README.md). Its [append-only reconciliation](../pmdi_hydroxyethyl_source_control_2026-10-26/registry_reconciliation.json) preserves the registry version used by this matrix and confirms that its four entries were unchanged.
