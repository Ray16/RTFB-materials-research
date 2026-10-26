# NDI external-control starting geometries

[Preparation summary](summary.json) · [Input plan](inputs/plan.json) ·
[Proposed state manifest](inputs/manifest.csv) · [Independent review](independent_review.json) ·
[Tests](tests_v2.xml) · [Experimental source](../ndi_external_challenge_2026-10-26/README.md)

The source-matched neutral NDI graph was sampled with ETKDGv3, seed 20261026,
followed by MMFF94 geometry optimization. Of 32 requested embeddings, 18 passed
the embedding-stage pruning and all 18 MMFF optimizations converged. Ten
post-optimization seeds were retained. Independent reconstruction checked the
declared connectivity and stereochemistry, finite coordinates, atom counts and
pairwise geometry separation; the smallest symmetry-aware heavy-atom RMSD is
0.176 Å. A second run reproduced the plan, manifest and all ten coordinate files
byte-for-byte under the recorded RDKit version.

The input layout uses `conformers/conf_*.xyz` and a separate three-state manifest
matching the project input schema. The manifest has **not** been merged into the
project library. An [isolated UMA adapter](../ndi_external_uma_2026-10-26_v2/README.md)
is now staged; a production dispatcher and DFT handoff remain pending. Charge and
multiplicity entries are proposed starting hints, not calculated state assignments.

These are neutral-graph, gas-phase force-field seeds. They are not charged-state
minima, quantum energies, solution free energies, a converged ensemble, redox
potentials or experimental validation. The source's quasi-reversible peak
midpoints remain provisional comparison quantities. A later calculation must
use the exact graph, state-specific optimization and native physical checks,
including alternative-spin assessment. This control is separate from the
post-functionalized polymer synthesis route.

The checked inputs are stored here for inspection; raw preparation and replay
locations are in [summary.json](summary.json). The earlier pre-layout draft is
preserved separately as [summary_v1.json](summary_v1.json).
