# Publication release preparation

The [updated local reviewer candidate](review_candidate_02.json) contains six checked
payloads totaling 818,081,211 bytes: committed scientific source, the current
manuscript upload and slides, historical review packages, electronic-control
checkpoints and the corrected molecular CPU replay. It is located at
`/nfs/lambda_stor_01/homes/rzhu/redox_release_candidates/2026-11-03_jcim_review_candidate_portable`.
Its manifest hashes all 18 bundled files and records the full validation build
of source commit `5ceb4e3b`. The earlier candidate and producer are preserved.
This is ready for local review; the
new full raw archive, public deposition and final metadata remain outstanding.

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
the 120.3-GB payload will fit. A suitable durable destination is pending user
input. No full raw archive has been started or marked complete.

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
