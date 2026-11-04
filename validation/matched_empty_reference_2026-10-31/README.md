# Matched empty-backbone references for the exact viologen panel

Two empty five-repeat references now have verified actual SCF starts under the
corrected implementation, one for each controlled backbone configuration. Their
exact graph identities match the empty members of the declared ionic occupancy
panel. All three reactive sites retain covalent benzyl chloride; sites 1 and 5
remain plain styrene. These are method-reference calculations, not a new
synthetically eligible redox candidate or an assumed experimental graft loading.

The warm-start selector checks every stored review/geometry hash, graph identity,
configuration, charge/spin, reactive-site pattern and solvent-model definition.
It selects the lowest stored 50-cm-1 G geometry within the declared old finite
batch. All twenty old states per configuration are considered. That old G is
used only to select coordinates. Neither old E nor old G enters a corrected
occupancy cycle, and the reference is freshly optimized and evaluated.

The seed evidence and selected coordinates are copied into each private
namespace. An initial attempt failed before QM because its source plan referred
to a current-workspace path hidden by the immutable archive bind. The failed
receipts are preserved. The namespace preflight now verifies every actual
namespace-visible source hash before importing the scientific runtime or taking
a GPU lease. The corrected attempts pass preflight and actual process/start,
CPU/thread, exclusive-residency and SCF-progress checks in [v2_startup.json](v2_startup.json).

The [summary](summary.json), seed lineages and [manifest](manifest.json) preserve
this scoped progress. The running optimizations, later Hessian minima and native
G remain unqualified; no balanced occupancy or conformer-ensemble conclusion is
claimed. The [78-test report](tests.xml) includes explicit neutral-reference
scope, changed-byte/graph/model rejection and the existing geometry, native and
comparison checks. Source hashes and as-used producer snapshots are retained.

Use `redox.workflows.matched_polymer_reference` with an explicit original
preparation, checked execution agreement, old finite-batch analysis, empty case
key and new output namespace. The standard probe/gated launcher performs
submission with the existing bounds. No missing property is filled from a warm
start or from the preceding numerical implementation.
