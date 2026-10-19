# Matched phenothiazine parent panel

Four exact free parents from Attanayake et al. (Chem. Mater. 2019,
DOI 10.1021/acs.chemmater.8b04770) are registered with three charge states each.
The [plan](plan.json) binds the structures, conformer seeds, charge/multiplicity,
scientific producers and protocol. [Source curation](../../reviews/2026-10-19/ptz_benchmark_sources/source_review.json)
provides five strict experimental target potentials and one strict spacing.
The other second-wave values remain excluded from accuracy scoring; their
calculated properties, if physically admissible, are predictions only.
B(MEEO)EPT is not scheduled because the paper's table and discussion disagree.

The [launch](launch.json) submitted 12 workers: eight V100s on lambda1 and four
on lambda2. [Actual start verification](start_verification.json) inspected live
controllers, one resident quantum process per GPU, disjoint two-CPU allocations,
inherited thread caps and node load. Startup does not prove completed properties.

Each worker starts from up to three FF-ranked ETKDG conformers (seed 0xC0FFEE),
uses charge-specific xTB-ALPB warm starts, and optimizes with the current full
grid-response DFT/SMD driver. UMA is not required for this supported multi-seed
path. It then computes/qualifies the active SP and native solution thermochemistry,
including the corrected internal-mode partition. All physical gates are unchanged.
Shared molecule locks exclude the existing exclusive repair/adoption paths;
state workers write disjoint directories. There is one attempt per state.
Failures preserve their outputs and require separate review before any repair.

Live status and logs are under `tasks/<parent>__<state>/`; numerical records are
under `calcs/dft/<parent>/<state>/`. Neither zero exit nor an output file admits
a value into the benchmark. The next step is independent physical/provenance and
energy-difference review, followed by a fresh table/benchmark assessment against
the matched Fc reference. No parent measurement is assigned to a grafted candidate.
No new accuracy result has been claimed at launch.
