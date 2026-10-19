# Matched-parent comparison after targeted spin qualification

Eight formal/half-wave potentials from five exact parents and four studies are
available under the stated state-specific spin models. The corrected native Fc
reference is 4.3791959016 V. RMSE is **0.18464 V**, MAE **0.16106 V**, and mean
calculated-minus-measured bias **+0.08431 V**. Model values use three-decimal
serialization, matching the earlier benchmark convention. No correction is fitted.

Both PMDI dianion comparisons are complete: the alternative triplets are higher
by 0.69162 and 0.70074 eV at the active SMD electronic level. These pass the stated
conservative electronic spin gate; their alternative thermal corrections are not
measured or asserted. The other states retain explicit DFT, UMA-large-gap or
standard-doublet confidence in `comparisons.csv` and `states.json`.

The first catholyte point, free EPT, is +0.537375 V calculated versus +0.278 V
measured. Its second wave is excluded because chemical reversibility is not
established in the primary source. Measurements belong to these exact free
parents and are not inherited by grafts. The supporting electrolyte is not
explicit in the continuum model; the PMDI source salt assignment is ambiguous.
This sparse, correlated set does not establish exact-graft accuracy or kinetics.

The independent raw-record reader reconstructs native free energies, reference
and potentials; 24 comparisons pass with stable input hashes, including measured
CSV inputs. No quantum SCF is repeated by this review. The new spin calculations
are separately recorded in `results/validation/benchmark_imide_spin_2026-10-19/`.
The older seven-point snapshot and superseded thermal/reference slide remain
historical. General benchmark and uncertainty consumers still need the same
completed-spin gate before refreshed global ranking/calibration tables publish.
