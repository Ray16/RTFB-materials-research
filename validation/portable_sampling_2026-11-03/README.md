# Portable manuscript validation

The complete article/SI validation and build now pass from a relocated source
tree. The strongest [filesystem-isolation check](home_isolation_02.json) rebases
the checkout to `/tmp/work`, hides the original user home except the existing
read-only Python runtime, disables networking and exposes no CUDA devices.
Only the shared CPU allocation lock is retained for resource coordination.
All numerical checks and compilation succeed: 21 article pages and 36 SI pages.
This is a workspace-independence test using the existing checked runtime, not
a fresh installation or a new quantum calculation.

This was subsequently repeated from the
[committed source export](source_export.json), not only the proposed overlay.
The [full build](committed_source_validation.json) and
[standalone sampling replay](committed_example_validation.json) both pass the
same home-hidden, network-disabled isolation for commit `5ceb4e3b`. The
[updated reviewer candidate](../publication_release_2026-11-03/review_candidate_02.json)
contains that exact source version and its validation record.

The [PDF comparison](rebased_pdf_equivalence.json) finds identical text and
rendered pixels on all 57 pages relative to the previously inspected PDFs.
No reported scientific value, figure or interpretation changes.

## Exact inputs and preserved checks

The [portable sampling example](../../examples/occupied_sampling_replay/README.md)
ships the bytes for all 183 direct analysis bindings, including twenty XYZ
geometries and the old producer sources. The reader requires the complete
original path/hash mapping and a matching analysis fingerprint. It verifies
each mapped input and refuses an unmapped historical path. Explicit current
producer files are recorded separately; they cannot replace the old producer
bytes. The shared sampling kernel still checks every state/context, source,
graph, chirality, duplicate and finite-partition condition.

The original ionic preparation checker is also bound by historical evidence
and remains byte-identical. Its maintained portable reader is
`redox.validation.ionic_preparation`. It permits only an explicitly declared
preparation-root relocation, verifies the original selected geometry strings
and hashes, and reproduces all twelve graph/XYZ/mass records and four balanced
cycles. The build separately verifies the old checker snapshot and binds the
new executing reader; only path/producer metadata may differ from the frozen
result.

[Thirty-six focused checks](focused_checks_02.xml) cover the shared analysis and
preparation gates, changed/missing input bytes, unsafe mappings, producer-version
separation, exact preparation relocation and rejected scientific mutations.
The [re-export replay](reexport_replay.json) reproduces all nine duplicate/rotor
scenarios from another fresh packet. The failed 2-kJ/mol sampling diagnostic,
eleven incomplete cells and lack of balanced occupied-ensemble admission are
unchanged.

## Evidence history

The proposed code/data overlay is recorded in `source_overlay_03.json` against
scientific commit `7d19efe8`; unrelated worktree edits are excluded. Earlier
attempts are preserved. They exposed the absolute preparation path and then
the immutable original-checker binding. The Python-audit build passes in
`isolated_build_03.json`; `filesystem_isolation_02.json` additionally masks
the original scientific directories for native libraries and child processes.
The final home-isolation run is stronger and also hides other original user
data. Earlier failed sandbox setups did not reach numerical validation.

The full 120.3-GB declared raw archive and public deposition remain separate
release work. This small portable input packet does not replace that archive
or claim to repeat SCF, Hessian or wavefunction validation.
