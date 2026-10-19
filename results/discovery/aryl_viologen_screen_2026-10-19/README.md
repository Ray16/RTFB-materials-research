# Diaryl-viologen computation

This campaign advances the nine [previously curated O-grafts](../aryl_viologen_2026-10-18/README.md)
into isolated gas-phase UMA screening. The declared comparisons retain both
N–aryl bonds: meta/para carboxyl substitution, dicarboxyl positional patterns and
two phenolic O-attachment sites. These extend the alkyl-cap series; they are
computational hypotheses, without inherited parent performance.

The [preparation](preparation.json) verifies nine unique inputs, 27 charge-state
definitions and 270 seed/state checks. Each graph has ten deterministic conformer
seeds. The exact isolated acid/phenol stoichiometry is held fixed across the
+2/+1/0 ladder. Equilibrium protonation in MeCN, exact monoattachment synthesis,
ester/ether durability and source-condition mapping remain unestablished.

The [UMA plan](uma_plan.json) freezes all inputs, model identity, code and 27
points. The 70-atom AV-3c graph supplies the three-charge resource pilot on reserved
16-GB T4s. Its [independent review](pilot_review.json) now passes all three charges,
with a 1473-MiB peak device allocation. Original pilot owners are terminal. The
remaining 24 points are assigned to 20 T4s and four H100s in the
[disjoint screen assignment](screen_assignment.json).
Every subsequent point retains measured memory, force/convergence, exact graph,
spin-parity and hydrogen-attachment gates. This qualification never admits DFT
or Hessians to 16-GB GPUs. All free suitable devices are used up to the declared
independent work after its dependency passes.

Per-point status is in `tasks/<id>__<state>/status.json`; full warm-start results
are in `uma/<id>/<state>/`. The gate launch manifests and host observations
distinguish submitted processes, actual GPU starts and verified records.
One actual attempt is allowed per frozen state; failures are preserved, and
three consecutive failures stop a worker. Resume requires inspecting original
process identities rather than clearing claims by age.

The additional [hydrogen-topology check](../../../src/redox/validation/protonation.py)
prevents a proton-transferred structure from silently passing this isolated
screen as the declared acid/phenol. Six regression cases and all 550 software
tests pass. The [stored-geometry audit](../../validation/hydrogen_topology_2026-10-19/README.md)
separately checks current records and reconciles warnings with existing integrity
exclusions. No earlier property was changed by this check.

All 27 gas UMA states have now passed the [terminal independent review](uma_completion.json).
All nine exact graphs are registered, with 27 library rows appended and the
verified conformer/UMA trees [installed byte-identically](registration_prepared/adoption.json).
The [solution campaign](../../validation/aryl_solution_2026-10-19/README.md)
has begun on suitable reserved V100/H100 GPUs. Its three-seed DFT mode uses FF
seeds and the UMA multiplicity; its actual seed selection is recorded separately.
New potential, usable-capacity and kinetic fields remain uncomputed. The source's
water/pH14/AgAgCl values are not MeCN/Fc labels for these exact grafts.
