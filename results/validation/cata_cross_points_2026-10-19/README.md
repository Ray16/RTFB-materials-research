# Cata-imide gas-relaxation diagnostics

Eight independent cross points ran on all eight free suitable H100s, with one
reserved worker per GPU, disjoint CPU pairs and bounded native pools. All eight
solver records completed; [actual starts](start_verification.json), [terminal
processes](terminal_verification.json), [source plan](plan.json) and
[review](review.json) separate execution from interpretation. These are gas SPs
at solution geometries; they do not validate solution kinetics.

| Exact graft | Event | Gas sum (eV) | QC exclusion |
| --- | --- | ---: | --- |
| Cata-NDI ethyl | First reduction | 0.4341 | Conformer change: 1.961 Å heavy-atom RMSD |
| Cata-NDI ethyl | Second reduction | 0.3547 | Gas-unbound anion; spin-ensemble kinetics unvalidated |
| Cata-ADI ethyl | First reduction | 0.3660 | Conformer change: 0.577 Å heavy-atom RMSD |
| Cata-ADI ethyl | Second reduction | 0.2846 | Gas-unbound anion |

The values remain diagnostics. All four retain flags and cannot serve as
QC-clean primary kinetic descriptors. Positive gas HOMOs do not invalidate the
separately admitted SMD reduction potentials. Conformer changes between
adiabatic solution minima also need a basin/ensemble treatment for reorganization.

The completed [full property-table review](../../reviews/2026-10-19/cata_complete_tables/review_status.json)
passes 32,492 independent arithmetic checks with stable inputs and zero
discrepancies. The full software suite passes 541 tests. Its new regression
ensures that per-record external evidence (such as saved qualification gradients)
is fingerprinted even when the native record JSON is unchanged.
