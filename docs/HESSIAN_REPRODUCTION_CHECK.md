# Prospective full-Hessian reproduction comparison

This numerical comparison is specified before any of the three corrected
119-atom full-Hessian controls has written its final matrix. It supplements
the ordinary physical checks; it does not replace them or establish chemical
accuracy, native chemical G, or ensemble convergence.

The [configuration](../config/polymer_hessian_reproduction.json) fixes separate
absolute tolerances for center energy (1e-8 Ha), Cartesian gradient (1e-6
Ha/Bohr), Hessian elements (1e-7 Ha/Bohr²), internal frequencies (0.01 cm⁻¹)
and nuclear terms at each rotor cutoff (1e-6 eV), together with a 1e-6 relative
Frobenius Hessian tolerance. These concern reproduction of the same numerical
model, not uncertainty relative to experiment. The nuclear tolerance is about
0.0001 kJ/mol, well below the manuscript's chemical effect sizes. The original
1% sampled derivative, stationarity, spin and minimum gates remain unchanged.

Every input must first have a successful independent full-Hessian reconstruction.
The comparator verifies that review's actual source hashes and requires the same
retained geometry, original state identity and numerical protocol. It does not
compare a corrected control with an uncorrected control under this acceptance
test. Distinct controls and all pairwise comparisons are required. Changed
backing, mismatched states or duplicate runs fail explicitly. Differences above
the limits are written as failed comparisons and return a nonzero exit status;
the tolerances are not adjusted after seeing the results.

After normal completion and separate physical review, run:

```bash
python scripts/validation/compare_screening_polymer_hessians.py \
  --review FRESH_INDEPENDENT_REVIEW.json \
  --review RESTART_INDEPENDENT_REVIEW.json --out NEW_COMPARISON.json
```

The same check can include the independently reviewed H100 control, or a future
raw-response continuation. A source-prefix audit, CPU test or completed response
block alone cannot serve as a completed full-Hessian input. Current results and
remaining limits belong in [the results index](../results/README.md).
