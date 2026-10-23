# Portable primary-result review package

[Current status and archive hash](current.json) identify a 32 MB prepared archive
on shared storage, containing 234 raw states and a source snapshot. It is not
yet a public DOI deposition or the complete auxiliary SI raw-state release.

The relocated standard-library replay reconstructs all 39 primary benchmark
endpoints, 67 reference grafting reactions, 20 base-assisted values, two
relative-basicity controls and the twelve-state finite polymer partition.
It also replays frozen manuscript claim checks. A separate NumPy/ASE CPU
check reconstructs 233 Cartesian-Hessian projections and nuclear partitions.
An audit hook rejects access to the original checkout/calculation paths.
Corrupt bytes, a coherently rehashed wrong thermal identity, an unbalanced
reaction and an escaping path are rejected.

The package includes original record bytes, frequencies, geometries, Hessians,
source code, explicit environment records and an internal README. Extraction
can be followed by:

```bash
python3 -I -S PACKAGE/code/scripts/release/reproduce_review_package.py PACKAGE
```

With NumPy/ASE and all six native thread caps set to one, use
`PACKAGE/code/scripts/release/review_package_hessians.py PACKAGE` for the
optional numerical reconstruction. This does not rerun SCF or repeat every
original stationarity, spin or finite-difference gate.

The [qualification records](../../validation/release_package_2026-10-23/qualified_v2/summary.json) and
[archive member hashes](../../validation/release_package_2026-10-23/archive_check.json)
are separate from the source package. Every one of the archive's 2,234 files
matches its source SHA-256. The package remains frozen; later evidence needs
a new package.

Recreate with `scripts/release/prepare_review_package.py --out NEW_DIRECTORY`
in the original checked scientific environment. The Git source commit and
release-tool overlay hashes are explicit in the package manifest. Use
`scripts/validation/review_release_package.py` to repeat relocation and
negative checks. The repository's MIT code license is included. Data-deposition
licensing and author metadata still require confirmation before deposition.
No accessible polymer capacity, new physical QC or experimental accuracy is
established by an offline reproduction.
