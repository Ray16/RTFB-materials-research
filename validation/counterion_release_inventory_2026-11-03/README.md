# Declared source inventory for the chloride-contact release

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

Ten bindings remain unresolved: five `thermal_policy.json` files and five
`hessian_diagnostics.json` files for historical DBU neutral/protonated/chloride,
PMDI precursor and local benzyl-chloride calculations. Their expected hashes
and referring records are listed in the resolution report. This establishes a
gap in the checked release inventory, not that the bytes have been lost or
that the independently reviewed chloride quantities fail physical QC. Further
archive/package recovery is required before claiming complete provenance
closure. Unresolved JSON branches were not traversed.

`inventory.py` produced the first pass; `resolve_known_content.py` performed
the exact-content lookup. Both preserve existing reports instead of overwriting
them. Reproduction requires access to the recorded shared-storage sources;
this inventory itself is not a portable numerical replay or public archive.
No SCF, Hessian or experimental validation was repeated by these scripts.
