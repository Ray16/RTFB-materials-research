# CN-PMDI spin ensemble diagnostic

Both the selected singlet and alternative triplet dianion have current active
energies and independently admitted corrected native thermochemistry. A CPU-only
calculation using the independent thermal reader gives a thermochemical gap of
**+0.01918 eV** (triplet minus singlet), **32.16% triplet population** at 298.15 K,
and an equilibrium free-energy correction of **−0.009968 eV** relative to the
selected singlet alone. Spin degeneracy is already included in each native free
energy and is not counted again.

These numbers are a source-bound diagnostic, not a released primary potential.
The second potential remains blank because the primary readers still need a
consistent equilibrium-ensemble implementation. A populated second spin surface
also prevents claiming validated single-surface kinetics. No new quantum jobs,
production records or acceptance thresholds were changed for this diagnostic.

See [summary.json](summary.json) for exact source hashes and [review.py](review.py)
for the independent-reader calculation.
