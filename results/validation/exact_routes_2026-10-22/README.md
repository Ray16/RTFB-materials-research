# Checked descriptor-to-route handoff

The [record](summary.json) reports an actual successor-stage execution on all 837
exact designs from the completed molecular descriptor campaign. Two CPU workers
on lambda2 used disjoint one-core allocations and single-thread native pools.
They completed eight and six chunks respectively. Every design remains on hold
because no exact-candidate route evidence has been supplied. No quantum jobs or
new candidate registrations occurred.

The new builder reconstructs the exact precursor/product and a unique mapped
heavy-atom bond from the same terminal dummies used by the molecule builder.
It checks the product, precursor, local chloride site, charge/electron-count
contract and balanced reaction. This represents a capped graft site on a preformed
polymer, not an experimentally demonstrated synthesis or a whole-polymer model.

The successor consumes verified upstream receipts and exact task identities.
Route evidence must bind the specific precursor/product and mapped attachment
bond. The existing route gate continues to enforce access, selectivity, conditions,
balance, core retention, source and actual-condition exclusions. Strong descriptor
values or evidence for another graph cannot open the full-redox gate.

Large route catalogs can be supplied as NDJSON and indexed on node-local disk;
mutable SQLite files are not placed on shared storage. Every recipe and upstream
task identity remains traceable. Frozen handler and builder sources are included
in [source_snapshot](source_snapshot/); raw inputs, journals and actual resource
records remain at `calcs/scale_validation/2026-10-21_exact_routes/`.

```bash
PYTHONPATH=src python -m redox.screening.route_stage prepare \
  --descriptors COMPLETED_DESCRIPTOR_CAMPAIGN --campaign NEW_ROUTE_CAMPAIGN
PYTHONPATH=src python -m redox.screening.route_stage worker \
  --campaign NEW_ROUTE_CAMPAIGN
PYTHONPATH=src python -m redox.screening.route_stage verify \
  --campaign NEW_ROUTE_CAMPAIGN
```

Run each CPU worker under the project's node/resource preflight with native
thread pools capped at one. The stage is parallel and resumable through the
existing chunk executor. Selected-subset conformer and quantum adapters remain
separate work; completion of this cheap stage does not establish redox accuracy.
