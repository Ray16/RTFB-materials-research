# Polymer physical-evidence replay after the refactor

[Checked summary](summary.json) · [Declared records](manifest_explicit_cwd.json) ·
[Sampling figure](../../results/figures/polymer_models/native_sampling_coverage.pdf) ·
[Figure provenance](figure.json) · [Visual review](visual_review.json)

All **60 retained, independently checked neutral-singlet states** reconstruct
from their archived Hessians, gradients, displaced SCF evaluations and nuclear
partitions. The reader verifies **12,386 original source/input bindings** through
one declared checkout-prefix relocation. Files and stored hashes are not rewritten;
missing archive files cannot fall back to the current checkout. Historical relative
paths require an explicit original-working-directory declaration. A final forced
hash pass checks every source again before publishing the completed replay.

| Exact graph/configuration group | Checked native states |
| --- | ---: |
| All-same empty backbone | 20/20 |
| Alternating empty backbone | 16/20 |
| All-same, site 2 occupied | 14/20 |
| Alternating, site 2 occupied | 10/20 |

The maximum frequency reconstruction difference is 6.97×10⁻¹¹ cm⁻¹; the maximum
nuclear-term difference over 25/50/100 cm⁻¹ rotor settings is 1.50×10⁻¹³ eV.
The maximum reconstructed sampled Hessian-column relative error is 3.92×10⁻⁵,
below the original 1% gate. Exact graph, charge/spin, solvent-model, gradient,
minimum-frequency, derivative, source and retained execution-resource checks
remain separate requirements. Seventeen focused tests pass, including rejection
of altered files, incorrect displacements, nonfinite forces/spin and conflicting
backbone-configuration aliases.

For the previously duplicate-reviewed all-same empty-backbone set, independent
log-sum arithmetic from the replayed G values reproduces the archived 10-to-20
changes: −0.85, −1.03 and −1.07 kJ/mol at 25, 50 and 100 cm⁻¹. The occupied
groups remain incomplete and their distinct-minimum/ensemble analysis is not
supplied here. In the all-same site-2 subset, the lowest checked late-index
structure is 3.84 kJ/mol below the best checked input among 1–10. Only seven of
those first ten have checked G, so this is **not** a formal complete-batch
convergence test. All panels use their own lowest checked G; absolute energies
from different chemical graphs are never compared.

The calculation model is the separately identified smooth-surface acetonitrile
model. Its G values cannot enter the primary-model PMDI cycle. The figure does
not infer experimental loading, association free energies, populations or
electronic coupling. The sixty records are checked states, not sixty established
distinct conformational minima. Broader occupied-graph sampling remains open.

To reproduce on the shared archive, set `PYTHONPATH=src` and cap OMP, MKL,
OpenBLAS, NumExpr and xTB threads to one. Run the following with the redox
interpreter (NumPy and ASE required; no GPU execution):

```bash
python scripts/validation/replay_archived_polymer_native.py \
  --archive /nfs/lambda_stor_01/homes/rzhu/archive/0_redox \
  --prior-root /nfs/lambda_stor_01/homes/rzhu/0_redox \
  --manifest validation/publication_native_replay_2026-10-28/manifest_explicit_cwd.json \
  --out results/validation/NEW_NATIVE_REPLAY/replay.json
```

The output must be new. Partial progress is explicitly unqualified until the
final source check passes. Plot with
`scripts/plotting/polymer_models/plot_native_sampling_coverage.py --replay REPORT
--out results/figures/polymer_models/native_sampling_coverage` and visually review
the rendered figure. Raw-record paths and hashes are retained in the large
report at the location recorded in [summary.json](summary.json), outside Git.

This does not rerun SCF, geometry optimization or DFT Hessians; no new quantum
jobs were launched. It does not qualify the current GPU production environment,
complete missing states, prove global conformer convergence or close the public
data-deposition requirement. The source archive is still required; this is not
yet a standalone public data package.
