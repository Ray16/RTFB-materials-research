# Streaming and parallel execution evidence

The [summary](summary.json) separates tested components from remaining full-workflow
scale qualification. The numerical controls below pass; scientific prediction
accuracy is not inferred from executor completion.

| Control | Checked result |
| --- | --- |
| Million-task packing | 1,000,000 unique synthetic inputs; 75.9 s packing; 68,792 KiB peak RSS |
| Distributed synthetic execution | Eight workers performed work on lambda1/lambda2; 1,000 chunks; 9.58 s worker span excluding preparation/startup |
| Independent synthetic output audit | [1,000,000 exact arithmetic checks](million_audit.json), zero discrepancies; input/result identity and hashes verified |
| Exact molecular pool | 1,593 recipes become 837 distinct tasks; all recipe lineage retained |
| Molecular descriptor equivalence | [837 serial/parallel comparisons](family_audit.json), zero discrepancies |
| Raw file count | 4,030 for the million-task control; no per-task control directory |

Four CPU workers were requested for the small molecular test; two completed all
chunks before the second node began useful work. Its 0.707-s execution span is
not a four-worker scaling or end-to-end speedup claim. The synthetic control did
perform work on all eight workers across both nodes. Controlled scaling curves
and per-worker peak memory remain future measurements.

The [execution guide](../../../docs/HIGH_THROUGHPUT.md) documents the streaming
API, bounded local dedup index, exact task/stage identity, atomic ownership,
checkpoint journals and explicit CPU resource allocation. Unit tests include a
real child-process crash, reuse of its four checkpointed results, and refusal to
reclaim a live owner. An interrupted-worker experiment on the actual shared
filesystem remains required before claiming full distributed recovery qualification.

Raw synthetic files are on shared storage at
`calcs/scale_validation/2026-10-21_million_control/`; the real molecular control is
at `calcs/scale_validation/2026-10-21_real_family/`. The audits record output hashes
and the [source snapshots](source_snapshots/million/plan.json) bind the exact
executed producer. Bulk inputs, journals and results are kept out of Git.
Historical outputs remain readable through their frozen source snapshots;
current workers refuse to execute a different producer under the same plan.

These controls evaluate orchestration and bookkeeping descriptors. They do not
constitute million-molecule redox prediction, attachment-route validation, DFT
throughput, polymer interaction accuracy or synthesis eligibility.
