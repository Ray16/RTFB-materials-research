# Refactor inventory and approval proposal

Status: approved by continuing the refactor, and implemented in this tree. Destination follows the user's explicit `0_redox_refactor` spelling. Source commit: e792621a0f140a627684ea411865d07d6e75e905. Source remains read-only. No push is planned. The module map and numerical comparison are in `MIGRATION.md`.

## Current inventory

`source_inventory.json` enumerates every Python library module, its SHA-256, module description, top-level definitions and imports. It is a source map, not a claim that every function has been independently scientifically validated.

| Existing package | Role and proposed disposition |
|---|---|
| build | Candidate/standalone/graft construction, geometry, manifests, open oligomer topology and contacts. Retain scientific logic in structures. |
| core | IO, hashing, metadata, protocol and runtime paths. Retain; separate path configuration from scientific identities. |
| properties | Potentials, references, Nelsen and solvent reorganization, capacity, ensembles, integrity and grafting. Retain named conventions and separate half/sum quantities. |
| qm | DFT/SP, UMA, xTB, solvent response, Hessians, thermochemistry and checkpoint adapters. Split by responsibility as below. |
| screening | Exact candidate registration, route review, evidence, feasibility and uncertainty-aware Pareto selection. Retain gates; no aggregate weighted score. |
| workflows | Campaign episodes mixed with resource management, queues, checkpoints, publication and reference installation. Extract shared stage orchestration and declarative campaign definitions; retain distinct admission policies. |
| validation | Scientific admission, provenance, reference/spin checks mixed with historical namespace and publication audits. Separate reusable QC from release validation and migration-only historical adapters. |
| reports | Benchmark, outer-sphere and workbook presentation. Consolidate with table/figure publication adapters. |
| root cli/current | CLI dispatch and compatibility entry points. Preserve command names and defaults; route to new services. |

Canonical does not mean newest filename. No module is proven dispensable merely by its date/name. Episode removal requires a caller/contract audit; unique numerical or rejection behavior survives as a shared policy and test.

## Data flow and public contracts

Exact precursor/product definitions and evidence → route proceed/hold/exclude → state/conformer and graft/oligomer structures → protocol-addressed UMA/DFT geometry and SP/Hessian records → source-bound thermal and QC admission → matched-reference potentials/reorganization/reaction and finite-occupancy quantities → qualified tables/figures and portable replay packets. Pending values remain absent. Completed checkpoints do not imply physical qualification.

Preserve `redox`, `redox-uma`, and `python -m redox`; commands: workflows, uma, dft, motifs, uma-check, status, sp, native-thermal, seed-reference, energy-audit, review, thermal-postprocess, build, build-validation, expand, search, source-audit, resources, launch, run-gpu. Maintain plan-versus-execute defaults. Script retention follows command callers, replay exporters and reusable campaign/validation/report entry points, not a recursive copy.

## Proposed target

```
src/redox/
  cli.py, __main__.py
  core/                 # IO, protocol, identities, paths
  structures/           # exact graphs, grafts, states, conformers, oligomers
  backends/             # dft, uma, xtb, solvent, hessian, response
  thermochemistry/      # nuclear partition, qrrho, projection, model definitions
  properties/           # redox, reorganization, reactions, ensembles, capacity
  screening/            # routes, evidence, admission, Pareto
  quality/              # numerical/physical/source/reference checks
  execution/            # reservations, resources, ownership, checkpoint transport
  workflows/            # stage graph, campaign specs, recovery, collection
  reporting/            # tables, figures, release/replay export
```

Dependencies flow from core to scientific domains to orchestration/reporting. Backends must not import campaign orchestration. Heavy optional libraries load only on paths that need them. Editable installation keeps runtime dependencies unresolved as today. Retain config and consumed curated/reference data, with a per-file dependency manifest. Seed-reference must remain independently testable.

## QM consolidation

- `thermal`: preserve the GFN2/ALPB diagnostic model separately from native DFT. Move partition helpers to thermochemistry and execution to the xTB backend.
- `native_thermal`: split native Hessian computation, QC, nuclear terms and orchestration. Preserve existing equations, reduction order, constants and definitions.
- `native_thermal_repair`: shared recovery workflow with archived-source promotion; not an alternative thermodynamic definition.
- `projected_thermal`, `qrrho`, `vibrations`: coherent thermochemistry layer, retaining projected/unprojected and historical model distinctions.
- `thermal_source`, `relocated_native`: provenance compatibility readers in quality. Historical producer hashes remain immutable. No broad hash bypass or re-certification by fiat.
- `hessian_memory`: execution budget/certificate policy. `hessian_resume`: additive CDS reconstruction. `hessian_restart`: orbital/tensor-frame restart. `saved_hessian`: complete-matrix reuse. Unify IO/identity checks in backend Hessian/checkpoint services; keep these recovery strategies explicit.
- `response_resume`: verified response reuse; `response_diagnostic`: instrumentation. Share frame/payload validation but preserve fresh residual checks and diagnostic-only status.
- SCF, occupied and sampling checkpoints: common durable IO/ownership mechanisms plus typed payload-specific admission. Do not erase different scientific checks to reduce file count.
- UMA current/policy/identity/nitro: one backend selection interface, retaining stock/Nitro qualification and parity contracts.

Source-bound numerical producer snapshots may need a narrowly scoped immutable compatibility area. Its exact contents must follow dependency and certificate tests; the live architecture must not become a facade over the entire old tree.

## Curated tests and exact oracle

Organize tests by structures, thermochemistry, properties, quality, execution, CLI and replays. Preserve every numerical or scientific rejection assertion, including those hidden in dated tests. Merge bookkeeping tests only after identifying their assertions. Keep topology, precursor/product identity, charge/spin/reference compatibility, invalid-source rejection, resource reservation and checkpoint durability coverage.

Required numerical families include energy conversion, potentials/reference_redox, reorg/solution_reorg/Nelsen, qrrho/native/projected/solution thermal, rigid-mode projection, contracts, UMA parity/identity, cached thermochemistry, both exact viologen packets, finite viologen and occupied sampling. GPU-dependent parity remains explicitly unexecuted unless existing CPU fixtures cover it; no new GPU work is authorized.

Capture old outputs and environment first. Compare structured scientific values exactly, including float representations where applicable; keep raw stdout hashes as additional evidence. Path-only differences must be enumerated, never blanket-normalized. Self-contained frozen examples alone cannot prove migrated-library equivalence: also compare the migrated kernels against their frozen implementations using the same inputs.

## Migration sequence after approval

1. Establish installable core, source/data manifests and oracle harness.
2. Port structures and route gates, then property math and thermochemistry.
3. Port backend interfaces and provenance/checkpoint readers with unchanged numerical kernels.
4. Consolidate execution and workflow specs; preserve reservations, CPU affinity/thread caps and recovery gates.
5. Curate reporting/exporters, scripts and tests; rewrite docs and complete old-to-new MIGRATION.md.

At each subsystem boundary run applicable contracts and all five cached replays. Stop on any numerical mismatch. Commit verified milestones on main as Ray Zhu, without coauthor trailers; do not push. No manuscript, slides, bulk runtime outputs, dated validation directories or browser binaries enter this tree.

Approval requested: target package boundaries, semantic (rather than filename-based) consolidation, and this test/migration sequence. Full source contract baseline and definitive per-module drop decisions remain work for the approved migration; no passing claim is made for unrun tests.
