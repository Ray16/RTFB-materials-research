# Screening implementation and polymer numerical assessment

The [current assessment](summary.json) records a reproduced GPU4PySCF 1.8.1
screening defect and an isolated implementation of the upstream correction.
The correction is **not yet adopted for production**. The chemical results
remain provisional for release while its effect on converged calculations is
assessed. This is an additional numerical requirement for M2/M3/M4/M5/M6/M7,
not evidence that every stored chemical result has a measured error.

The official 39-atom water-cluster control uses an all-ones density matrix,
def2-TZVPD, revPBE and sorted/unsorted grids. It is a fixed-density integration
test, not a molecular reaction or an SCF solution. The source and correction are
preserved with their upstream URLs and Apache 2.0 license in
`tests/fixtures/screen_index_v1_8_1.*` and `tests/fixtures/licenses/`.
The installed wheel's binary agrees with its distribution RECORD; its SM70
disassembly retains the original one-barrier loop. The run-local correction
adds the upstream read/write barrier without editing the installed library.
[Upstream issue](https://github.com/pyscf/gpu4pyscf/issues/919),
[correction](https://github.com/pyscf/gpu4pyscf/pull/925).

| Device | Installed XC grid-order difference (Ha) | Corrected difference (Ha) |
| --- | ---: | ---: |
| V100 | −1.018303459e-3 | +1.055923349e-9 |
| H100 | −1.018303410e-3 | +1.055923349e-9 |

The corrected density integrals also pass the original absolute 2e-9 bound.
Independent reviews reconstruct the differences and verify closed ownership,
exclusive residency, source hashes and resource bounds. The separate synthetic
mask fixture agrees with its independent CPU bound for all sixty comparisons
per device; it did not reproduce the original defect. See
[V100 review](V100_screening_review.json), [H100 review](H100_screening_review.json)
and their copied integration inputs. A prior NVRTC compilation failure and a
RAM prelaunch refusal are retained outside the checkout; neither is a failed
physical control.

## Chemical derivative controls and pending impact

The [closed-shell OH-anion control](V100_closed_review.json) uses r2SCAN-D4,
def2-SVPD and SMD MeCN at a fixed geometry. All six Hessian columns at both
0.0005/0.001 bohr steps pass the unchanged 1% relative gate, for installed and
corrected implementations. The largest column error is 1.2645e-5. The largest
installed/corrected differences are 1.42e-14 Ha in energy, 8.31e-15 Ha/bohr in
gradient and 6.62e-12 Ha/bohr² in Hessian. This small system does not establish
the correction's effect on polymers or benchmark redox energies.

The OH-radical control failed SCF convergence at a displaced geometry and is
preserved as a failed control. The distinct
[methyl-radical control](V100_open_review.json) passes all twelve Cartesian
columns at both steps for both variants. Its maximum relative Hessian error is
3.7883e-4; installed/corrected energy, gradient and Hessian differences are
2.84e-14 Ha, 1.28e-14 Ha/bohr and 1.74e-12 Ha/bohr². These two small controls
qualify only their named tests. H100 derivative controls have experienced slow GPU-driver
queries; actual process identities and calculation records determine progress.
Four independently owned V100 polymer controls completed two installed and
two corrected repeats on the exact three saved geometries around the previously
failed atom-31/z column. Their [startup receipt](polymer_start.json) verifies
disjoint CPU/GPU assignments, source bindings and thread caps. Four separate
closed reviews and the [implementation comparison](polymer_comparison.json)
verify the final force components. The largest repeated/variant force difference
is 1.16e-9 Ha/bohr; the earlier approximately 2e-6 difference is not reproduced
in these trials. CDS forces are identical. Comparisons against the original
Hessian column pass, but a full corrected polymer Hessian and the defect's
chemical-cycle impact remain unqualified. Raw source-frozen namespaces are under
`/nfs/lambda_stor_01/homes/rzhu/redox_kernel_controls/2026-10-31_polymer_*`.

The two old dispatch controllers were stopped at exact child-free boundaries;
their [native](native_admission_hold.json) and [electronic](energy_admission_hold.json)
receipts preserve PID/start ownership. No scientific GPU workers were signaled.
New uncorrected production admissions remain held. Installed packages, the
original scientific source tree and numerical acceptance thresholds are unchanged.

## Separate native repeat and response recovery

A full native repeat passes the existing protocol after the original 1.2138%
column failure. [Comparison](comparison.json) and [independent physical replay](native_replay.json)
retain both attempts. Fifteen identical SCF coordinates have energy differences
at most 7.28e-12 Ha and gradient differences up to 1.99e-6 Ha/bohr. Their Hessians
differ by at most 1.67e-8 Ha/bohr². The failed coordinate was atom 31/z;
N29–C38 identified force components, not the displaced coordinate. The repeated
maximum column error is 0.7044%, below the original 1% gate. Other columns retain
step-dependent residuals. The cause and finite-difference step convergence are
not established by passing that gate.

The [census](native_census.json) contains 122/240 physically reconstructed input
states under their stored protocol, including this repeat. Complete occupied
ensembles and balanced ensemble cycles remain missing. These counts do not
qualify the new implementation or update the article's separate 75-state cohort.

The full checkpointed V100 response execution control has
[passed independent file/frame/resource review](response_review.json): seven
blocks cover 315 Cartesian right-hand sides, with final mo1/e1 arrays saved.
Peak GPU memory was 19,050 MiB and peak host RSS 16.49 GB. Source-bound solver
records meet their original residual targets. The review does not independently
reevaluate those GPU residuals. Original grid arrays are absent, so this is not
a qualified exact Hessian restart and supplies no native G. It used the installed
screening implementation; it cannot qualify the correction.

## Reproduce

Use the locked scientific interpreter, `PYTHONPATH=src`, one CPU for preparation
or review and all native thread pools capped at one. Set
`NUMPY_MADVISE_HUGEPAGE=0` before importing NumPy. The versioned entry point
`python -m redox.workflows.screening_control` provides source-frozen prepare,
node-specific probe and guarded launch stages. `--water-reference` selects the
upstream integration control; `--molecular-control closed|open` selects full
Cartesian derivative controls. `prepare-polymer --spec ORIGINAL_NAMESPACE
--screening-variant installed|corrected` retains the existing three-point
scientific diagnostic in a new namespace.

GPU launches require the lambda-fleet workflow, exact dependency/resource
preflight and `gpu_reserve`. These bounded controls have separate documented
output budgets; the production queue's 64-GiB storage floor remains unchanged.
Closed controls are reviewed with `review_screening_control.py`,
`review_screening_molecular.py` and `review_response_checkpoint.py` in
`scripts/validation/`. They refuse changed inputs and never infer thermochemistry
from execution completion. `manuscript/check_screening_controls.py` reconstructs
the frozen integration arithmetic with standard Python and no CUDA.

M1, M4, complete M7 and final M9 remain open. The screening implementation's
chemical impact is now an explicit additional release requirement.
