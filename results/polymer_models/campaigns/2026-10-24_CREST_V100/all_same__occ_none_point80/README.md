# Second backbone configuration: bounded V100 resource point

The all-same empty-backbone CREST search has finished, and twenty exact starts
are staged separately from the twenty alternating-configuration starts. This
branch evaluates only the first all-same conformer at fixed coordinates.
It computes gas/SMD energies and a geometry-level gradient as a resource test;
it establishes no optimized minimum or native free energy.

The screened full/range-separated DF arrays require about 75.8 GB before the
declared conservative host margin. This separate experiment enforces an 80 GiB
RSS ceiling and 32 GiB free-RAM floor. [plan.json](plan.json) retains the input,
source, graph and physical settings. [launch.json](launch.json) and
[startup_verification.json](startup_verification.json) record actual execution;
[resource_point/status.json](resource_point/status.json) owns completion.
Raw records are under
`calcs/validation/2026-10-24_CREST_V100/all_same__occ_none_point80/`.

Broader DFT re-ranking, matched ensemble G and sampling convergence remain pending.
