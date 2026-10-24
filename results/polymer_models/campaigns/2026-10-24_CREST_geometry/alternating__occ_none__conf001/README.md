# First final-CREST geometry refinement on V100

This branch optimizes the first of twenty staged conformers from the completed
alternating-configuration empty-backbone search. The search supplied 4,329
topology/stereochemistry-checked structures. Search hit counts are not population
weights. The model has five repeats, three unreacted benzyl chloride handles and
two styrene padding units; it contains no grafted motifs.

An [independently reviewed exact-input resource point](../../../../reviews/2026-10-24/crest_empty_V100_point/summary.json)
precedes this separate geometry stage. [plan.json](plan.json) binds the exact
input, graph, protocol and control evidence. [launch.json](launch.json) and
[startup_verification.json](startup_verification.json) record the reserved V100,
disjoint CPUs, native limits and actual scientific start. The
[geometry status](geometry/status.json) owns completion. Host RSS is guarded at
80 GiB; GPU peak, graph/stereochemistry, optimizer convergence, stationarity and
independent total-energy derivatives still have to pass.

Coordinates, SCF and optimizer checkpoints are durable under
`calcs/polymer_models/2026-10-24_CREST_geometry/alternating__occ_none__conf001/`.
Canonical hybrid energies and native thermochemistry follow separately. One
running geometry does not establish DFT re-ranking of twenty starts, a broader
ensemble free energy or sampling convergence.
