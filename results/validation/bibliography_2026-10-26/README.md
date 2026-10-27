# Current manuscript bibliography and context audit

[Final metadata audit](checked_final/summary.json) ·
[Five non-DOI primary sources](non_DOI_primary_review.json) ·
[Primary transport-source assessment](transport_source_review.json)

The current bibliography has 48 entries. All 43 DOI-bearing references match
registered title, first-author and publication-year metadata. The five non-DOI
entries were directly checked against the primary patent records, the Tartu
collection and the official software sites. The US patent inventor's name was
corrected to Souâd Boulmaâz. Exact versions and source revisions remain separate
from software website publication dates.

The original parallel registry pass hit HTTP 429 limits. Its incomplete results
are preserved; they were not promoted. The final producer uses one network
request at a time, bounded Retry-After handling, immutable metadata-cache hashes
and source snapshots. It refuses a bibliography that changes during its run.
A false author mismatch caused by a TeX accent command was corrected in the
normalizer and covered by a regression check. Ten bibliography checks and six evidence-capture checks pass.

RedDB supplies context for systematic aqueous redox enumeration, QCArchive for
quantum-data infrastructure, and primary Oyaizu and Akhoury studies for useful
charge propagation in redox polymers. Their factual registered metadata accompanies the portable checked audit. Raw registry replies and primary article abstracts remain in shared storage rather than the versioned release.
The current paper keeps aqueous and nonaqueous scope distinct and does not
transfer transport rates or mechanisms from these polymers to our neutral
finite models. RDKit is now cited at the graph-processing step.

DOI identity checks do not establish the correctness of every content claim,
complete author lists or page ranges. A registry record for the def2 paper gives
only its first page; the existing full page range is retained. The current
editorial changes do not validate new candidate properties, Nitro parity,
interaction energies or polymer ensembles. The previous reference-status file
is preserved separately; `manuscript/references_status.json` points to this audit.

The final pass additionally cites all seven experimental benchmark studies. Six primary studies were previously present in the endpoint provenance but absent from the bibliography. A separate [source inventory](../benchmark_source_inventory_2026-10-26/current.json) accounts for every unchanged endpoint. Registry author names were normalized to Unicode NFC for pdfLaTeX; the original and rendered bibliography remain canonically identical.
