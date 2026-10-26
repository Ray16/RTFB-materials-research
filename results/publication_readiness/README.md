# Publication preparation

The pipeline has substantial checked evidence, but publication readiness remains
incomplete. The [requirements and current evidence](current.json) are a work plan,
not a release certificate. The [readiness guide](../../docs/PUBLICATION_READINESS.md)
defines what must be demonstrated before each requirement can close.

The [consolidated local reviewer bundle](../releases/manuscript_bundle/README.md)
now combines the source/evidence snapshot and four overlapping raw-state
packages, with relocated standard-library and locked NumPy replays. It is not
a public archival deposit; data licensing, author metadata and DOI remain open.

| Evidence | Location |
| --- | --- |
| Fresh software checks | [Test records](software_checks.json) |
| Observed compute environment | [Package versions](compute_environment.json); not a portable lock |
| Fixed retrospective experimental panel | [39 endpoints, 30 molecules](../validation/publication_benchmark/panel_v1.json) |
| Independent benchmark coverage and errors | [Current analysis](../validation/publication_benchmark/current.json) |
| Completed native benchmark campaign | [Frozen completion plan](../validation/publication_benchmark/native_completion_2026-10-21/plan.json) |
| Verified benchmark worker starts | [Resource and process checks](../validation/publication_benchmark/native_completion_2026-10-21/start_verification.json) |
| Polymer physical/cycle validation | [Independent report](../polymer_models/workflow_validation.json) |

The benchmark calculations are exact-parent method controls. Measured parent
properties do not establish accuracy for new grafts or whole polymer beads.
Raw native outputs remain in `calcs/dft/<parent>/<state>/` on shared storage;
completed molecule receipts bind their actual output hashes in the completion
campaign's `done/` directory. Failed attempts and partial outputs are preserved.

Regenerate a fresh independent benchmark review after calculations finish, then
revalidate its evidence before using the numbers in a manuscript. Current
reports are scoped snapshots and can lag ongoing calculations.

The fixed retrospective panel has 39/39 admitted endpoints and 111 independent arithmetic checks with no discrepancies. The [22 October publication checkpoint](../reviews/2026-10-22/benchmark_complete_publication/summary.json) preserves its earlier 14-page draft. The [current manuscript](../../manuscript/README.md) includes the finite PMDI cycle and the later DIPEA UV/NMR source checks. Broader grafting/interaction validation and the complete public quantum release remain open.
