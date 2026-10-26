# Complete finite PMDI cycle: portable raw-state addendum

[Archive and hashes](current.json) · [Extracted replay with original reads blocked](relocation.json)
· [Adverse-input test results](negative_tests.xml) · [Replay instructions](../../../examples/neighbor_cycle_replay/README.md)

The 3.95 MB compressed package includes fourteen original cached states, adding
the completed paired PMDI chain to the unchanged first-site package. It contains
the four matched chain geometries and their electronic records, native Hessians,
modes and nuclear partitions, plus shared local, reagent and redox/Fc records.
All 94 archived files pass the extracted inventory/hash check. With the original
checkout inaccessible, isolated standard Python reproduces fourteen balanced
reaction comparisons and the complete conditional cycle at 25/50/100 cm^-1,
including both attachment orders. No scientific or production package is loaded.

Seven integration checks pass, including the nominal replay and rejection of
scope promotion, NaN G, duplicate state identity, a changed/rehashed paired energy,
a nonfinite/rehashed cutoff result and a path escape. The tests accept an explicit
extracted package through `REDOX_NEIGHBOR_REPLAY_PACKAGE`; they skip when the raw
archive is absent instead of depending silently on the original cluster paths.

The arithmetic replay does not reintegrate nuclear terms or repeat the original
physical QC. It supplies no conformer convergence, isolated motif binding,
electronic coupling or synthesis admission. The original first-site archive
remains unchanged. This addendum complements the v5 and numerical-control
packages; it does not merge their overlapping states into a new unique count.

The archive and extracted data remain on shared storage at the paths in
`current.json`. Public data licensing, final author metadata and DOI deposition
remain open. This is a checked local reviewer artifact, not a public deposition.
