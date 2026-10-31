# Benchmark implementation sensitivity and full-Hessian controls

All **71 electronic states underlying the 39 benchmark endpoints** now pass
independent review after fresh calculations with installed and corrected
screening. The same exact geometries, charge/spin assignments,
ωB97M-V/def2-TZVPD/SMD MeCN protocol and original multi-guess/Fc occupation
strategy are used. Each selected energy is backed by its actual SCF checkpoint;
normal ownership, residency, runtime and resource checks pass.

| Comparison | Largest state-energy difference (Ha) | Largest electronic contribution to a relative-potential shift (V) |
| --- | ---: | ---: |
| Installed versus stored | 7.340e-10 | 1.666e-8 |
| Corrected versus stored | 7.344e-10 | 1.667e-8 |
| Corrected versus installed | 2.274e-12 | 8.662e-11 |

The [standard-library replay](energy_replay.json) reconstructs every difference
and Fc-referenced electronic shift from the [raw scalar inputs](energy_inputs.json).
The correction therefore does not explain the benchmark bias through its
electronic-energy contribution in this fixed-geometry panel. Nuclear terms
have not been recomputed for all benchmark states, and no new qualified
potential is assembled from mixed implementation components.

The largest RKS/UKS resource pilots have 1,350 AOs. Their measured GPU peaks
are 27,920/28,280 MiB and host RSS about 1.45 GB, below their prospective
32-GB GPU/40-GiB host bounds. These independent checks release the remaining
panel through a resumable, source-bound dispatcher. It discovers every
suitable free policy GPU, checks CPU/RAM commitments, retains uncertain starts
and scientific failures, and applies a three-failure node circuit breaker.
Prelaunch refusals are separate from scientific attempts. The full panel's
declared durable-output bound is 6,338,051,216 bytes, not a measured final size.

Both [closed-shell](H100_closed_review.json) and
[open-shell](H100_open_review.json) MeCN derivative controls also pass on H100,
complementing the V100 controls. These are fixed OH-anion and methyl-radical
method controls, not minima or material decisions. All Cartesian columns at
0.0005/0.001 bohr pass the original 1% threshold for both implementations.
The largest relative column error across these four device/state controls is
3.789e-4.
The [molecular raw inputs](molecular_inputs.json) and
[standard-library replay](molecular_replay.json) reproduce all 144 saved
Cartesian derivative checks without NumPy or quantum software.

## Full polymer Hessian and checkpoint scope

The initial 119-atom V100 control saved its SCF, partial solute Hessian, solute
H1 and original grid arrays before a proactive response-workspace refusal.
The installed aligned batching formula requires about 3.05 GiB for its
1,193-AO minimum block; the earlier 3-GiB ceiling cannot meet that minimum.
The revised preparer computes the required whole-GiB ceiling from actual
AO/occupied dimensions before launch. Its 4-GiB cap remains bounded by actual
available memory. Physical thresholds and the production storage floor are
unchanged.

One unqualified restart prototype incorrectly replaced the outer solution H1
method with the saved solute component. It was stopped after exact ownership
verification, and its [rejection receipt](rejected_restart.json) is preserved.
No Hessian or free energy from that prototype is accepted. Corrected reuse
intercepts the solute method underneath the PCM/SMD wrapper, which recomputes
the solvent H1 term. A regression test requires that term to survive. A separate
fresh calculation supplies the full-matrix comparison; successful quantum
restart qualification remains pending. Future controls also checkpoint the
complete solution H1 explicitly, with its component scope and source frame.
The [saved-component check](PCM_H1_component_check.json) finds a nonzero PCM
contribution, with a largest H1 element of 0.1145 Ha/bohr. This confirms that
the corrected restart retains the missing term; complete physical Hessian and
restart qualification still require the fresh-run comparison and final QC.

The raw controls are outside the checkout:

- `redox_kernel_controls/2026-10-31_full_Hessian_V100_fresh4/`: fresh full control.
- `redox_kernel_controls/2026-10-31_full_Hessian_V100_resume_PCM4/`: corrected component reuse.
- `redox_kernel_controls/2026-10-31_full_Hessian_H100/`: initial architecture control.

These paths are under `/nfs/lambda_stor_01/homes/rzhu/`. Their exact process
identities, guarded starts and live stages determine execution; source files
or a sleeping dispatcher alone do not prove progress.

The earlier seven-block/315-coordinate response execution probe has a narrower
scope than a native response calculation: its saved H1 is solute-only and its
operator uses frozen-solvent response. It remains an instrumentation and
execution diagnostic. Its array/resource checks are valid within that scope;
it does not qualify a full equilibrium-solvent response, an exact native
Hessian restart or native G. Original production native drivers call the
solvent-wrapped `Hessian.kernel()`, preserving both terms; the scope mistake
concerned the separate diagnostic/restart prototypes.

## Reproduction and remaining work

`scripts/fleet/run_screening_energy_panel.py` coordinates the fixed panel.
`src/redox/workflows/screening_energy_panel.py` prepares exact state inputs and
prospective size-based checkpoint budgets. The independent review and collector
are `scripts/validation/review_screening_energy.py` and
`collect_screening_energy_panel.py`. Electronic sensitivity replay uses
`src/redox/validation/screening_sensitivity.py` without quantum software.

Use `PYTHONPATH=src`, capped native thread pools and
`NUMPY_MADVISE_HUGEPAGE=0`. GPU execution additionally requires the lambda-fleet
preflight/reservation workflow. The scoped raw panel is at
`/nfs/lambda_stor_01/homes/rzhu/redox_kernel_controls/2026-10-31_benchmark_energy_panel/`;
small independent state reviews are versioned here. Bulk orbitals/checkpoints
remain on shared storage with hashes and locations.

Complete polymer Hessian/restart checks, balanced occupied ensembles, exact
candidate polymer assessment and final public deposition remain open. The
stored-protocol working census remains 122/240 inputs. This milestone does not
constitute complete reviewer resolution or submission readiness.
