# Bounded V100 point for the alternating empty backbone

The original twenty-start CREST queue and numerical protocol are retained in
the [plan](plan.json). This separate experiment evaluates only its selected
starting geometry. Both gas/SMD SCFs and the geometry-level gradient finish
within the enforced 80 GiB RSS/32 GiB free-RAM limits. The observed GPU peak is
20,390 MiB. The [independent resource review](../../../../reviews/2026-10-24/crest_empty_V100_point/summary.json)
checks finite energies individually, input geometry, normal pilot status, source
hashes, completed ownership and host limits.

The conservative storage forecast remains in the plan; it was not presented as
an observed allocation. This bounded test qualifies its named energy/gradient
resource point. Full DFT re-ranking, stationary minima, native thermochemistry
and the matched ensemble remain separate stages.
