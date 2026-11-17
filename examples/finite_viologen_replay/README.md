# Finite viologen cycles and motif comparison

Run with Python 3.10 or newer, from any working directory:

```sh
python -I /path/to/finite_viologen_replay/replay.py
```

No installed project, scientific packages, network, GPU or original calculation
paths are needed. The package retains byte-identical inputs from the paper's
3 November editorial evidence snapshot and verifies their hashes before use.

The replay reconstructs twelve cached electronic/nuclear component sums, all
four finite ionic occupancy cycles at three rotor cutoffs, and twelve frozen
motif triplets in gas and MeCN. It checks charge/atom-count balance, exact cycle
members, recorded state/QC, compatible model/thermal/reference conventions and
the motif XYZ hashes and A/B/AB concatenation. The existing independent motif
checker is retained verbatim; its sole project import, file hashing, is replaced
with the equivalent standard-library adapter. No chemical calculation is changed.

The main-cutoff G cycles are approximately -3.94, +22.58, +8.28 and -5.30 kJ/mol.
Frozen methyl-capped motif electronic interactions are positive in all four
frames, while whole-chain electronic nonadditivity is negative in two. The
output keeps those different quantities explicit.

This replays cached arithmetic and recorded context. It does not rerun SCF,
wavefunction, Hessian, mode or experimental validation, nor establish an
equilibrium ensemble, binding G, coupling, actual loading or material admission.
Original raw paths remain provenance strings only. Full raw-data deposition is
separate from this small replay package. For molecular redox and local
quaternization results, see the neighboring `exact_viologen_corrected_replay`
example.
