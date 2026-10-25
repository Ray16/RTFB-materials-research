# Matched-model hybrid energies after geometry review

[Forty declared entries](batch.json) cover the two empty-backbone configurations
and twenty starts per configuration. One identical checked reference is linked
directly; the other 39 entries require fresh hybrid gas/SMD energies. Deferred
preparation waits for each exact geometry's independent review. The measured
same-constitution/state/basis/model RAM class supplies only an execution bound.
No energy or geometry is copied to another conformer.

The [reference review](../../../reviews/2026-10-25/smooth_energy_point80_followup/checked/2026-10-25_smooth_energy_point80__all_same__occ_none__conf001.json)
records the actual 80 GiB host-RSS class and 32-GB-class GPU execution. Workers
retain durable gas/SMD SCF checkpoints and distinct experimental energy records.
Every completed case needs its separate independent energy/resource review.
The controller uses available suitable GPUs and the shared CPU/RAM-promise and
reservation checks. CPU dependency probes defer GPU imports until reservation.

Dynamic progress is in `dispatch.json`; initial and later scientific starts are
separate from accepted energy records. Raw calculations remain under
`calcs/polymer_models/2026-10-25_smooth_empty_electronic_batch/` and the original
reference case. All new plans and the batch bind their sources and exact inputs.
These electronic records do not supply native G or conformer-ensemble convergence.
