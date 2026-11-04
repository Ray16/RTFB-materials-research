# Working polymer collection: 115 reconstructed inputs

Fourteen additional phthalimide-control states pass an independent reconstruction
of their saved native Hessians, sampled derivative columns and nuclear free-energy
terms. The [replay](../../results/validation/polymer_collection_2026-10-31/replay.json)
verifies 34,299 source bindings. Its largest sampled-column relative error is
3.382e-5 (0.003382%); its lowest internal frequency is +7.43 cm^-1. All fourteen
retain the same explicitly identified smooth-CDS MeCN model. This is physical
record replay, not a new quantum calculation or an external-accuracy benchmark.

The [working census](native_census.json) now contains **115 of 240 declared
input states**. The
[coverage figure](../../results/figures/polymer_models/2026-10-31_collection/native_occupancy_coverage.pdf)
distinguishes the twelve exact graph/configuration groups:

| Configuration | Empty | Site 2 | Site 3 | Site 4 | Sites 2,3 | Sites 2,4 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| All same | 20/20 | 19/20 | 12/20 | 11/20 | 6/20 | 9/20 |
| Alternating | 20/20 | 17/20 | 1/20 | 0/20 | 0/20 | 0/20 |

Input counts include duplicate minima. Neither these counts nor individually
checked G values establish complete Boltzmann ensembles, balanced occupied-ensemble
cycles or an experimental occupancy distribution. Reviewer M4 remains open.
The article retains its separate, source-frozen 75-input numerical comparison;
an incomplete working collection does not replace that scientific analysis.

## Remaining calculations and storage

The retained raw outputs of the fourteen completed cases occupy approximately
0.679 GiB per single-graft state and 0.982--0.986 GiB per double-graft state.
Their largest files are first-derivative response arrays, not the final small
Hessians. The reusable
[measurement script](../../scripts/validation/measure_native_storage.py) records
the actual file inventories, namespaces, source hashes and observed sizes in
[storage_measurement.json](../../results/validation/polymer_collection_2026-10-31/storage_measurement.json).
Sixty single-graft and 65 double-graft native inputs remain uncollected. Applying
the observed class maxima estimates a 104.84-GiB full retained footprint for
those unfinished inputs. Some files are already allocated by running jobs, so
this is not an additional-write budget. It also does not bound temporary growth
or include upstream electronic outputs and other users' writes.

The volume had approximately 61 GiB free at measurement, below the unchanged
64-GiB admission floor. Reopening the gate alone would not provide enough space
to certify storage capacity for the cohort. A further writable shared destination
would allow the remaining campaign to proceed. No scientific record, frozen evidence or old checkpoint
was deleted or moved, and no new GPU calculation was launched in this collection.

Four older site-2 calculations are still confirmed live by exact PID/start
identity. Three short CPU/GPU samples show CPU activity and no sampled GPU
utilization. Their saved partial-Hessian/first-derivative matrices remain available,
but a completed orbital-response receipt or final Hessian has not appeared.
These observations do not identify the hot function or prove scientific
progress. The jobs were not stopped or restarted based on age. The
[liveness record](../../results/validation/polymer_collection_2026-10-31/old_site2_liveness.json)
preserves the observations separately from physical qualification.

## Reproduce

Run the redox interpreter with native thread pools capped at one. The CPU tools
allocate one disjoint core, preserve old outputs and verify source identities.
Raw scientific inputs remain on shared storage under
`/nfs/lambda_stor_01/homes/rzhu/redox_polymer_continuations/`. Their original
namespaces and hashes are explicit; access to those records is needed for full
physical replay. The released cached CPU examples are a separate portability
demonstration.

```bash
python scripts/validation/replay_namespace_native.py \
  --manifest results/validation/polymer_collection_2026-10-31/selection.json \
  --out /path/to/new/replay.json

python scripts/plotting/polymer_models/plot_occupancy_coverage.py \
  --census validation/polymer_collection_2026-10-31/native_census.json \
  --out /path/to/new/coverage

python scripts/validation/measure_native_storage.py \
  --census validation/polymer_collection_2026-10-31/native_census.json \
  --sample-replay results/validation/polymer_collection_2026-10-31/replay.json \
  --volume /nfs/lambda_stor_01/homes/rzhu/0_redox \
  --out /path/to/new/measurement.json
```

The census's `source_replays` list supplies all eight inputs to the maintained
`collect_native_census.py --replay ...` command. Its original 20-input-per-group
target and all physical thresholds are unchanged. The previous 101-input
[snapshot](../polymer_recovery_2026-10-31/native_census.json) is preserved.
