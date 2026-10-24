# Portable scientific-result review package

[Current status and archive hash](current.json) identify an approximately 32 MiB
prepared archive on shared storage with 243 thermochemical states, 26 accepted
electronic oligomer models and six frozen-fragment single points. The earlier
primary package remains preserved. No public DOI deposition is claimed.

The relocated standard-library replay reconstructs all 39 benchmark endpoints,
67 reference grafting reactions, 20 base-assisted values, two relative-basicity
controls, the twelve-state finite polymer partition, five neutral-basis
sensitivities and two restart energy comparisons. It also reconstructs three
interaction/deformation partitions from raw electronic records and checks frozen
fragment coordinates against their archived complex. Electronic-only records
remain explicitly without free energies; the failed stacked trimer is excluded.

An optional NumPy/ASE CPU check reconstructs 242 Cartesian-Hessian mode sets and
nuclear partitions, plus both restart geometry/frequency comparisons. Original
checkout/calculation opens are blocked. Seven negative cases reject corrupt
bytes, wrong thermal identity, unbalanced reactions, external paths, promotion
of an alternative-basis certificate, wrong electronic spin and promotion of an
electronic energy to a free energy. These are numerical/cached-record checks,
not fresh SCF or repeated original physical qualification.

The archive includes raw geometries, frequencies, Hessians, energy records, source
code, explicit environments and an internal README. Run after extraction:

```bash
python3 -I -S PACKAGE/code/scripts/release/reproduce_review_package.py PACKAGE
```

With NumPy/ASE and all six native thread caps set to one, use
`PACKAGE/code/scripts/release/review_package_hessians.py PACKAGE` for the
optional numerical check. [Qualification records](../../validation/release_package_2026-10-23/qualified_v5/summary.json)
and [all 2,498 archive-member hashes](../../validation/release_package_2026-10-23/archive_check_scientific_v5.json)
are separate from the frozen package.

Recreate with `scripts/release/prepare_review_package.py --out NEW_DIRECTORY`
in the original scientific environment. A Git commit plus exact producer/test
overlay hashes defines the code. `scripts/validation/review_release_package.py`
repeats relocation and adverse-input checks.

Execution/numerical diagnostic artifact closure, fresh GPU restoration, data
deposition licensing, author metadata and public deposition remain open. The
repository's MIT code license is included. No new material performance, physical
admission or experimental accuracy follows from this review package.

A [GitHub draft release](draft_release.json) now holds the checked v5 archive,
package manifest and checksum file. All three assets were downloaded again and
verified byte-for-byte. The draft requires repository-maintainer access and is
not a public raw-data deposition or DOI. Its numerical epoch remains
`2026-10-23_all_axes_checked`; current polymer work and subsequent manuscript
changes are outside that snapshot. Data-reuse terms, author metadata and the
final publication snapshot remain unresolved.

A [separate fresh compatible quantum deployment](../../validation/gpu_environment_restore_2026-10-23/README.md) is now independently checked for two complete spin controls. It is outside the v5 archive and does not alter that frozen package, create a DOI or establish general platform portability.
