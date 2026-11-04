# Publication release preparation

The [complete local reviewer package](review_candidate_06.json) is finalized
inside the existing canonical release directory:
`results/redox_release_candidates/2026-11-03_jcim_review_candidate_submission`.
Use its `CURRENT.md` and `manifest_5f50c61_with_raw.json`. All earlier manifests,
packages and source versions remain unchanged.

The seven payloads include the shared 15-page Article and 41-page SI, aligned
submission files, self-contained manuscript source, scientific source commit
`6629fcd1`, existing scoped replay/control and slide packages, and the complete
verified declared raw archive. The final manifest binds 37 package files.
The scientific source archive matches all 12,377 Git blobs.

The [isolated source check](current_source_validation.json) reconstructs the
scientific-evidence documents and compiles the shared manuscript from that
committed archive with original user data hidden and networking disabled.
Its shared PDFs' page text matches the checked 15/41-page version. Compiler logs
are included. Extraction occurred in private tmpfs with no persistent work folder.
The installed runtime was reused; no new quantum calculation was performed.

The [full raw verifier receipt](full_raw_archive_verified.json) confirms every
archive member and the final compressed checksum. The 107,289,881,732-byte archive
contains 46,789 unique research-evidence objects (120,381,992,561 uncompressed
bytes), plus five manifest/support members. The source bytes were hashed during
creation and the entire archive was independently read back. This establishes
byte integrity for the declared inventory, not new physical qualification or
coverage of implicit dependencies. Runtime dependencies and third-party
publications remain explicitly external.

Local archive preparation and verification are complete. Public deposition,
data reuse metadata and final author declarations remain pending. No journal
submission or public data deposit is claimed.

## Source replay and resolved portability gap

The [source export](source_export.json) is exactly commit `7d19efe8` and excludes
unrelated uncommitted work. From its 11,414 extracted files, the
[isolated checks](source_replay.json) replay the finite viologen cycles/motifs,
corrected molecular cohort and phenothiazine oxidation with original scientific
paths blocked. These checks reconstruct cached numerical results, not fresh QM.

The initial [whole-manuscript validation attempt](portability_gap.json) identified a
dependency on absolute source paths in the occupied-sampling reader.
The first attempt stopped at a nested CPU-allocation guard in the test harness;
the second used the builder's own allocation and reached the actual relocation
failure. Neither is counted as a successful full validation build. The separate
self-contained TeX upload still compiles to the checked 21/36-page PDFs.
The [subsequent portable-input fix](../portable_sampling_2026-11-03/README.md)
preserves those source and physical-context checks. A fresh committed-source
export now passes the complete validation build with original user data hidden,
networking disabled and the source rebased to `/tmp/work`, using the existing
checked runtime. Its standalone sampling replay passes the same isolation.
The updated review candidate includes this source and its matching validation
record; the failed initial attempts remain historical evidence.

## Raw inventory and storage

The [seed manifest](raw_seed_manifest.json) selects frozen publication evidence,
including the 31-state sampling snapshot rather than a live campaign census.
Recursive named hash maps also reach historical calculations and restart
records. Inclusion does not qualify a failed or unfinished state.

The [first inventory](raw_inventory_summary.json) verified 45,464 files and
identified 7,991 unresolved bindings. Two explicit namespace-resolution passes
([first](namespace_resolution_summary.json),
[second](namespace_resolution_summary_02.json)) locate all declared bindings:
53,640 physical files and **120,269,364,790 bytes of unique content**. Recorded
namespace mounts are used only as location declarations; no mounts or jobs are
executed. Every newly located file must match its original SHA-256. Two worker
bindings resolve to the unchanged control script; the old analyzer is recovered
byte-for-byte from Git. The first resolver source is preserved separately.
Prior verified raw files are reused by the resolution passes rather than being
described as freshly rehashed. Final archiving must recheck all included bytes.

The complete inventories remain on shared storage at the paths and hashes in
the summaries. Named `source_sha256`, `evidence_sha256` and `producer_sha256`
maps are covered; implicit runtime dependencies and other schemas are not.

The [storage observation](storage_headroom.json) found 14 live calculation owners,
28.16 GB of remaining output commitments and a 32-GiB reserve, leaving about
51.36 GB for new archive writes at that observation. The
[compression pilot](compression_pilot.json) measures four full H1-family files;
the two large arrays retain 95.5% and 96.1% of their size at zstd level 3.
This sample is not a full-archive bound and provides no basis for assuming that
the 120.3-GB payload will fit. At that observation, a suitable durable destination was pending user
input and no full raw archive had started. Storage has since become available;
see the current archive section below.

## Maintained packaging tools

`redox.validation.source_inventory` supplies shared declared-map traversal;
`scripts/release/stage_declared_sources.py` stages a resolved object package,
and `scripts/release/stream_source_archive.py` can stream it without a duplicate
uncompressed tree. Both use the same inclusion and relocation rules. Runtime
dependencies and third-party publications stay explicit and external. The
standalone `verify_source_objects.py` checks extracted objects and support code;
the streaming verifier checks every archived member, rejects unsafe/duplicate
members, and requires the complete expected member set.

[Sixteen checks](archive_checks.xml) cover exact relocation, multiple source
versions, missing or wrong content, altered support code, corrupted resumed
objects, source changes during archiving and changed manifests. Archives
preserve the 32-GiB floor and use one zstd worker. The dated scripts in this
directory are audit producers, not calculation-submission entry points.

## Current raw addendum and release blocker

The [addendum](raw_addendum_summary.json) extends the resolved baseline to the
current seven-contact evidence epoch: 54,101 physical files and 120,434,748,131
bytes of unique content, with no unresolved declared bindings. Fourteen explicit
seeds cover the 13 changed snapshot contents and the seven-outcome cohort. The
shared inventory reader hashes all 461 newly discovered files. The prior 53,640
file observations are explicitly inherited, not described as current rehashes.
The [checkout producer check](inherited_producer_observation.json) separately
rehashes 326 current-checkout text/code dependencies and finds no changed bytes.
Implicit dependencies and live campaign outputs outside the frozen seeds remain
outside this inventory's scope.

[Nineteen tests](incremental_inventory_checks_02.xml) verify inherited-versus-new
versions, wrong copies, late relocation, incomplete baselines, altered aliases,
object staging and streaming archive integrity. Archive creation still rechecks
every included byte; incremental inventory is not an archive-integrity certificate.
Use `scripts/release/extend_source_inventory.py` with the hashed baseline summary
and explicit seed manifest for future additions. The original reader is preserved
as `source_inventory_v1.py`; earlier inventories and bundles are unchanged.

The [fresh storage observation](storage_headroom_02.json) finds 78.45 GB free.
The included research objects require 120.38 GB uncompressed, plus a 34.36-GB
(32-GiB) reserve: a minimum 76.29-GB shortfall before active calculation promises.
The existing compression pilot does not justify assuming the full archive fits.
That observation preceded the newly available storage. A full raw archive
is now running as described below. Author details, reuse licensing, a permanent release identifier
and public deposition remain outstanding. This small inventory addendum accompanies
the immutable reviewer bundle through the results index; it does not rewrite that
bundle or its scientific-source commit.

## Full raw archive completion

The [current archive inventory](current_raw_archive_inventory.json) resolves the
canonical source locations and preserves the previous inventory. The one changed
producer path was recovered byte-for-byte from Git and stored as
`raw_builder_source.py`; its original SHA-256 remains unchanged. The archive covers
46,789 unique research-evidence objects totaling 120,381,992,561 bytes, with runtime
dependencies and third-party publications kept explicitly external.

`full_raw_archive_owner.json` records the live process identity and archive path.
The writer uses one zstd worker, hashes every source while streaming, and preserves
the 32-GiB storage floor. It then independently verifies the complete archive
member set and hashes. `full_raw_archive_verified.json` now records a passing result for all 46,794
members and the compressed archive checksum. The gated finalizer admitted the
archive only after that result and owner-stage closure.
The [13 stream/package checks](raw_archive_stream_checks.xml) include existing
metadata directories and refusal of conflicting filenames. No new release or
working directory was created.

The final manifest gate requires an exact matching inventory and manifest hash,
complete member/object counts, the verifier's final archive checksum, and an
unchanged local file size. Pending or mismatched receipts cannot admit the raw
archive into the final package. The original version-five manifest remains intact.
