# Corrected primary exact-viologen CPU replay

Run `python replay.py replay --package .` with the hash-locked NumPy dependency.
The ten-state cohort contains nine newly repeated primary molecular Hessians,
all nine saved-gradient panels and analytic monatomic chloride. It reconstructs
two Fc-referenced potentials, three precursor-specific quaternizations and
disproportionation at all three rotor cutoffs. The runner reuses the original
frozen example implementation, changing only its required FD-panel count to nine.

Only cached numerical data are read. No SCF/GPU/xTB calculation, original project
path or installed redox package is required. Actual SCF orbitals, resource review
and independent spin reconstructions have separate source-bound receipts/raw
storage; they are not rerun by this NumPy cache. Code, modes, partition equations,
geometry and energy/protocol records are included. No experimental/speciation,
conformer-ensemble or material admission or full historical executable closure
is inferred. The original six-panel example remains preserved separately.
