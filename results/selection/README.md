# Candidate decisions: all required axes first

Use [current.json](current.json) for the current publication, its immutable input
snapshot and machine-readable decision/front records. The
[selection specification](../../docs/ALL_AXES_SELECTION.md) defines the rule.

The separate [PMDI linker matrix](../discovery/pmdi_linker_matrix_2026-10-26/README.md) applies the same conjunctive policy to four route-held graph hypotheses. Its two new products are not included in the frozen 879-graph ledger.

The separate [hydroxyethyl PMDI source control](../discovery/pmdi_hydroxyethyl_source_control_2026-10-26/README.md) applies the same policy to a proposed O-graft with a competing imide N–H; its route remains held.

The separate [hydroxyethyl phenothiazine proposal](../discovery/ptz_hydroxyethyl_source_2026-10-26/README.md) applies the same policy to a held exact graft outside the frozen ledger.

The later [TEMPO exact-case assessment](../grafting/tempo_exact_support_2026-10-26/README.md)
uses the same conjunctive policy but is not part of the frozen 879-graph ledger.
Its route passes and its overall material decision remains pending. The source
registry used by the historical ledger is preserved in the linked reconciliation.

| Output | Scope |
| --- | --- |
| [Descriptor Pareto figure](../figures/selection/descriptor_pareto.pdf) | 56 qualified-for-display molecular models; eight members of the joint E1/E2/nominal-capacity descriptor set |
| [All-axis evidence matrix](../figures/selection/candidate_axis_matrix.pdf) | Pending/pass/fail decisions; supporting model/proxy information is separate from a criterion pass |
| [Acceptance policy](../../config/candidate_feasibility.json) | Exact identity, model scope, uncertainty/sensitivity and context requirements |
| [Focused expansion plan](../../docs/TARGETED_CANDIDATE_EXPANSION.md) | 34 design slots across seven cores, plus an attachment control; not 34 generated/qualified molecules |

The fresh October 23 candidate snapshot contains 75 candidates and one control.
Its independent reconstruction has 628 checks and zero arithmetic discrepancies.
The union with the systematic pool is 879 distinct graphs: 878 candidate designs
and one comparison control. All 878 candidate decisions are pending in this
profile. No full-axis shortlist has been published. Several evidence quantities
and four application targets remain unresolved; this is not evidence that all
878 chemistries are poor.

The descriptor plot is a projection of a three-objective set, so its outlined
points are not joined into a purported two-dimensional boundary. No candidate
uncertainty interval is imputed. Molecular nominal capacity is not presented as
functionalized-polymer capacity. Its current inputs and front-member IDs are in
[descriptor_pareto.json](../figures/selection/descriptor_pareto.json).

## Reproduce or refresh

```bash
export PYTHONPATH="$PWD/src"
export OMP_NUM_THREADS=1 OMP_THREAD_LIMIT=1 MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 XTB_THREADS=1
python scripts/pipeline/publish_candidate_decisions.py --out NEW_PROJECT_RESULT_DIRECTORY
```

The publisher checks the fixed 39-endpoint benchmark, refreshes the independently
audited candidate feed if its source inventory is stale, publishes all-axis
records and both figures, verifies source stability and only then updates the
current pointer. It makes no new QM submissions. Both the normal finalization
workflow and candidate-overlay update call this successor stage. Existing
molecular scorecard/Pareto outputs are retained as descriptor analyses rather
than silently reinterpreted as a full candidate shortlist.

Raw state records remain under `calcs/dft`, `calcs/grafting` and
`calcs/polymer_models`. The current decision publication stores its exact input
and producer snapshots, source hashes and publication manifest. Older frozen
manuscript evidence remains unchanged.
