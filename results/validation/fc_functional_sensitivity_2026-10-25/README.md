# Fixed-geometry Fc functional diagnostic

[Independent review](independent_review.json) reconstructs six completed gas/SMD
electronic calculations: neutral and oxidized Fc with omegaB97M-V, B3LYP and
PBE0, each using def2-TZVPD and the same qualified solution geometries. The
original-functional SMD energies reproduce the primary records within 1e-6 eV.
All six SCF/spin, runtime and exclusive-residency checks pass.

Relative to omegaB97M-V, the SMD electronic gap decreases by 0.11840 eV with
B3LYP and 0.06581 eV with PBE0. Holding the primary target/nuclear terms fixed
increases the benchmark bias from +0.18959 to +0.30799/+0.25540 V. These
reference-only counterfactuals are diagnostics, not alternative validated redox
protocols. No geometry, nuclear partition or alternative potential is computed;
the primary reference and benchmark remain unchanged.

[Plans](campaign.json) bind the input and producer identities. Raw SCF records
are under `calcs/validation/fc_functional_sensitivity_2026-10-25/`, one directory
per functional/state. The independent review binds their hashes. The immutable
paper snapshot `manuscript/evidence/2026-10-25_reference_diagnostic/` includes all
six electronic JSON records and reconstructs the gaps and counterfactual metrics.
Broader method accuracy and unique attribution of the reference error remain open.
