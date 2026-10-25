# Smooth-model hybrid-energy and RAM pilot

[The source-bound plan](all_same__occ_none__conf001/plan.json) scores the exact
independently checked empty-chain geometry with fresh omegaB97M-V/def2-TZVPD
gas and MeCN SMD single points under the distinct smooth-CDS model. It creates
`experimental_energy.json`; it cannot supply a primary energy record, native G
or candidate admission. Gas and solvated SCFs have durable checkpoints.

The preceding attempt at
`results/polymer_models/campaigns/2026-10-25_smooth_energy_pilot/`
exceeded its 32 GiB host-RSS guard during density-fitting allocation. Its own
guard stopped that process; the attempt and its source remain intact, with no
accepted energy. That geometry-stage RAM bound was insufficient for the hybrid
stage. The [existing same-constitution hybrid resource point](../../../reviews/2026-10-24/crest_empty_V100_point/summary.json)
measured about 74.8 GiB and passed an 80 GiB bound. This separate restart uses
that bound with 32 GiB free headroom and retains the same electronic equations,
coordinates, GPU ceiling and physical criteria. The legacy point supplies only
resource evidence; its geometry/energies are not transferred to this model.

[Current checked status](../../../reviews/2026-10-25/smooth_followup_checkpoint/summary.json)
and the [dedicated independent reviewer](../../../reviews/2026-10-25/smooth_energy_point80_followup/status.json)
distinguish placement, actual execution and accepted outputs. Raw records remain
under `calcs/polymer_models/2026-10-25_smooth_energy_point80/`, with failed-attempt
and measured-resource hashes bound into the plan. A successful energy/resource
pilot would still leave native thermochemistry and ensemble convergence open.
