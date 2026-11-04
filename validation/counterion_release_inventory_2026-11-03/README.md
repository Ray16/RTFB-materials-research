# Declared source inventory for the chloride-contact release

## Resolved inventory and local package

All ten remaining metadata bindings now resolve to exact copies under
`archive/0_redox/calcs/grafting/`; the
[mapping and hash checks](grafting_metadata_resolution.json) preserve their old
logical names. The metadata introduces no additional dependencies in the three
declared hash-map fields. All 163 initially unresolved bindings are therefore
located. This completes the declared-map inventory, not an inventory of every
implicit runtime dependency.

The [local package verification](package_verification.json) checks 5,528 unique
objects (3,720,557,167 bytes), supplying 7,207 file mappings and the included
relocations. Every object passes an isolated standard-library verification with
original scientific paths blocked. The package is at
`/nfs/lambda_stor_01/homes/rzhu/redox_release_candidates/2026-11-03_declared_contact_sources`.
Its manifest records the 19 runtime dependencies and one third-party publication
as external; they are not redistributed. This is a local raw-evidence package,
not a fresh numerical replay, public deposit or license grant.

The maintained staging and verification commands are
`scripts/release/stage_declared_sources.py` and
`scripts/release/verify_source_objects.py`. Staging retains a 32-GiB free-storage
margin and can explicitly resume an unfinished objects-only directory by
rehashing its contents. Initial staging correctly stopped when an external
dependency appeared in the relocation list; included and external relocation
maps are now separate, and the preserved partial objects passed rechecking.
[Six focused tests](release_checks_02.xml) cover relocation without original
paths, corruption, unsupported mappings, unsafe object names, symlinks and
the external-dependency/resume boundary. `check_package.py` records the actual
full-package isolated check.

## Original inventory and resolution history

The [inventory](inventory.json) starts with the six checked contact outcomes,
their association ledgers and matched primary graft/chloride parents. It
recursively follows top-level `source_sha256`, `evidence_sha256` and
`producer_sha256` path-to-hash maps in verified JSON. Aggregate provenance also
reaches historical PMDI, base and other molecular records, so the resulting
payload is broader than six contact calculations. Implicit dependencies and
other metadata schemas are not claimed to be covered.

The pass verifies 7,217 files: 7,070 unversioned research/source files, 127
versioned project artifacts, 19 runtime dependencies and one third-party
publication. Their total size is 3,772,480,614 bytes; content deduplication within
those technical categories gives 3,743,152,667 bytes. These categories do not
assign ownership or redistribution rights. Bulk copying and deposition have
not occurred.

An explicit checkout-to-archive mapping resolves 3,774 path/hash bindings.
The first pass leaves 163 absent logical bindings. The
[content-resolution pass](content_resolution.json) locates and rehashes exact
copies for 153 of them among already verified files, retaining both original
paths and replacements. Different versions at the same logical path remain
distinct hash bindings; existing certificates are never rewritten.

Ten bindings remained unresolved at that stage: five `thermal_policy.json` files and five
`hessian_diagnostics.json` files for historical DBU neutral/protonated/chloride,
PMDI precursor and local benzyl-chloride calculations. Their expected hashes
and referring records are listed in the resolution report. The subsequent
grafting-path recovery above closes that gap without changing the independently
reviewed quantities. These JSON files contain no further declared hash-map
branches.

`inventory.py` produced the first pass; `resolve_known_content.py` performed
the exact-content lookup. Both preserve existing reports instead of overwriting
them. Reproduction requires access to the recorded shared-storage sources;
this inventory itself is not a portable numerical replay or public archive.
No SCF, Hessian or experimental validation was repeated by these scripts.
