# Consolidated manuscript reviewer bundle

[Archive and status](current.json) · [Standard-library replay](stdlib_replay.json) ·
[Neighbor Hessian isolation](neighbor_native_isolation.json) ·
[Numerical-control isolation](controls_native_isolation.json) ·
[Preparation specification](../../../config/manuscript_reviewer_bundle.json)

The local 44.9 MB archive gathers a committed source/evidence snapshot and four
previously checked raw-state packages: the v5 primary/SI reviewer package, the
first-site addendum, the numerical-control addendum and the complete finite
PMDI neighbor-cycle addendum. Their molecular identities overlap; state counts
must not be added. Third-party source PDFs/XML are not included or silently
relicensed; the manuscript records their DOIs and scope.

The bundle was extracted outside the original checkout. With original project
reads blocked, isolated standard Python reconstructed the primary benchmark
and reactions, first-site sums, electronic controls, complete neighbor cycle
and current manuscript claims. A separate locked NumPy 2.2.6 environment
reconstructed 13 neighbor Hessian mode sets and the monatomic partition, and
eight numerical-control geometry/Hessian records. The source tar and all nested
archives passed their file inventories and SHA-256 checks. This does not rerun
SCF, stationary-point, spin or original experimental validation.

To reproduce the small CPU example after obtaining the archive, extract the
outer bundle and the neighbor addendum with `zstd` and `tar`, then run:

```bash
python3 -I -S /path/to/2026-10-26_neighbor_cycle_native_replay/reproduce.py \
  --data /path/to/2026-10-26_neighbor_cycle_native_replay
```

The output includes an exact PMDI potential, fourteen balanced comparisons and
the conditional neighbor cycle. Each nested package contains its own README;
the neighbor package also describes its hash-locked NumPy Hessian replay.
The producer is `scripts/release/prepare_manuscript_bundle.py`, and the
relocated standard-library verifier is
`scripts/release/validate_manuscript_bundle.py`. Bulk archives and the relocated
replay remain on shared storage at the paths in [current.json](current.json).

This is a **local reviewer artifact**, not a public deposition. The data reuse
license, author metadata and archival DOI remain unconfirmed. A source tar and
cached replay do not establish broad GPU-environment portability, complete
conformer sampling, experimental accuracy or candidate admission.
