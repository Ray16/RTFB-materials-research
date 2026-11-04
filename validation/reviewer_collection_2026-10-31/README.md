# Consolidated local reviewer collection

A checksum-verified **422,312,575-byte** collection contains the committed
`66254a7b` source/evidence snapshot and eight previously checked raw-data/replay
packages. It includes the primary scientific package, first-site and neighboring
cycles, numerical controls, TEMPO/local-motif addenda, the 75-state numerical
cache and exact-viologen CPU example. Record identities overlap; their counts
must not be added as unique states.

The archive is
`/nfs/lambda_stor_01/homes/rzhu/redox_release_candidates/2026-10-31_reviewer_collection.tar.zst`.
Its SHA-256 is
`dd7ea1b9ee9279020cb5e48941b807cd392f997b9424d1bdf3eeb91b6ed4c1c2`.
The [checked summary](summary.json) and [collection manifest](collection_manifest.json)
retain source-commit, nested archive locations, hashes and explicit scopes.

Assembly verifies all eight original archive hashes and tests every compressed
stream. The outer archive is checked member-by-member against the collection
manifest. A fresh extraction of the committed numerical replay then reconstructs
71 electronic states, 39 reference-cancelled electronic shifts and all 144
molecular derivative checks using isolated standard Python, with original
project/archive/raw-calculation reads denied. This repeats the new cached
numerical scope; it does not claim that every nested package's scientific
replay was repeated during collection assembly. Their earlier checked archives
remain byte-identical.

The source and eight packages remain separately extractable, with their own
locked environments and replay instructions. The collection includes completed
committed evidence; current full-Hessian and charged-optimization pilots are
not promoted into its accepted scientific records.

Author/affiliation metadata and a data-license choice are still needed for DOI
deposition. No public DOI, final deposition or submission readiness is claimed.
Third-party experimental publications retain their original rights. The
collection does not change the open polymer, exact-candidate or final narrative
requirements.

The versioned assembly and check entry points are
`scripts/release/reviewer_collection.py` and
`scripts/release/check_reviewer_collection.py`. They use bounded CPU allocation,
streaming compression, preserved source packages and a new output destination.
