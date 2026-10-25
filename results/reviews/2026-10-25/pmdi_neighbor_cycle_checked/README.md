# Exact PMDI conditional grafting cycle

[Independent reconstruction](summary.json) closes the four-state cycle for the
same neutral five-repeat backbone: empty, site 2, site 3 and sites 2/3. All four
states have admitted electronic energies and independently reconstructed native
nuclear partitions. Atom/charge balances and both attachment orders close.

At the primary 50 cm-1 cutoff, nonadditivity is **−1.8616 kJ/mol**. Defined
DBU/separated-ion grafting at site 3 changes from −96.9582 kJ/mol on the empty
chain to −98.8198 kJ/mol with site 2 grafted. The reverse order gives the same
direct double-grafting sum. The base/precursor/byproduct contribution cancels
from nonadditivity. Cutoffs 25/50/100 give −1.0698/−1.8616/−1.8208 kJ/mol;
these settings are not statistical uncertainty bounds.

The paired minimum has core centroids 16.30 Angstrom apart and a 69.0-degree
interplane angle. Thus this is a finite whole-chain grafting response for
separated cores, not a measurement of close-contact stacking strength, isolated
core binding, electronic coupling or an electron-transfer rate. Each graph
contributes one minimum; broader conformer/configuration convergence and material
eligibility remain open. Absolute grafting values retain the DBU/speciation
qualifications of the earlier first-site analysis.

[Conditional-energy figure](../../../figures/polymer_models/pmdi_neighbor_grafting.pdf),
[actual geometry](../../../figures/polymer_models/pmdi_neighbor_geometry.pdf)
and their adjacent JSON caches bind the source data. The paired raw state is
under `calcs/polymer_models/2026-10-25_residual_native/`; the summary binds
all four raw state locations and hashes. The [completed component review](../pmdi_pair_native_residual_cuda_class_followup/component/summary.json)
reports a +5.77 cm-1 minimum mode, 7.62e-6 Ha/Bohr maximum gradient and
3.43e-6 maximum relative sampled-column error. The earlier local publication
handoff error is preserved separately; the checked component is now registered
with the [matched baseline](../../../validation/pmdi_occupancy_2026-10-23/baseline_current.json).
