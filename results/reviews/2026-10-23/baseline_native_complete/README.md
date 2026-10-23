# Independent complete two-pose baseline

All twelve native free energies reconstruct independently from their active
electronic energies, saved DFT Hessians, projected internal frequencies and
nuclear/spin partitions. All six pairs are structurally distinct under the
declared 0.05/0.10/0.20 Å duplicate criteria.

At the primary 50 cm⁻¹ rotor cutoff, adjacent and one-spacer nonadditivities are
−5.169 and −11.393 kJ/mol. The 25/50/100 cm⁻¹ sensitivity ranges are
[−5.169, −4.639] and [−12.337, −11.393] kJ/mol. These ranges are model sensitivity,
not confidence intervals. Conditional DBU-assisted grafting values are stored
separately for separated ions and contact pairs using the frozen reagent models.

The finite-set configurational lowering of each graph is at most 0.14 kJ/mol.
This only describes these two known poses. It does not bound the error from
unfound lower-energy conformers, establish a global minimum or justify a general
spacing rule. The wider searches remain necessary.

[Source-bound numerical review](summary.json), [producer snapshot](producer_source.py)
and [baseline figure](../../../figures/polymer_models/conditional_ensemble_baseline.pdf).
