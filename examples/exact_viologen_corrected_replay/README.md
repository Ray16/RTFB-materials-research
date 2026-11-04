# Corrected exact-viologen CPU example

Run `python examples/exact_viologen_corrected_replay/replay.py replay` with NumPy.
The hash-locked CPython 3.11/Linux requirement is in `data/requirements.lock.txt`.

Ten cached states reconstruct the exact molecular redox/quaternization cohort:
nine newly repeated primary-model Hessians, all nine saved-gradient panels,
and analytic chloride translation. All six declared cycles retain exact method,
state and reference definitions and three rotor settings. No GPU, SCF, xTB
library, original project installation or original calculation path is needed.

This is cached mode, derivative and thermochemical arithmetic replay. Original
SCF orbitals/resource checks and independent spin reconstructions have separate
source-bound receipts in the validation record; the NumPy example does not
repeat them or establish experimental accuracy, sampling convergence, speciation
or material admission. The earlier six-panel example remains intact at
`examples/exact_viologen_replay/`. Full raw-data deposition and final metadata
remain outstanding. The maintained exporter is
`scripts/release/primary_molecular_example.py`; the standalone runner preserves
the original numerical implementation with the FD-panel requirement raised to nine.
