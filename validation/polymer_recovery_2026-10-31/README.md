# Occupied-polymer reconstruction and recovery

The working census contains **101 independently reconstructed native input
states**, up from the article's frozen 75-state comparison. The
[census](native_census.json) and
[coverage figure](../../results/figures/polymer_models/native_occupancy_coverage.pdf)
distinguish all six occupancy graphs and both declared backbone configurations.
The site-2 sets contain 19/20 and 17/20 checked inputs. Other occupied sets are
incomplete. Input counts include duplicate minima; they do not establish a
Boltzmann ensemble, global sampling, a balanced occupied-ensemble cycle or an
experimental grafting distribution.

Twenty-four new occupied states and two previously completed site-2 states pass
additional Hessian/derivative/nuclear reconstruction. The three replay reports
verify 58,599, 1,631 and 2,594 source bindings. These are neutral phthalimide
controls at the same explicit smooth-CDS MeCN model; they are kept separate
from primary-model reaction components and exact-candidate material decisions.

## Scheduling recovery

An old foreign-residency refusal occurred before the scientific host guard or
any SCF started. Its leftover residency record kept the dispatcher treating
lambda13 as unsettled. Recovery now requires the exact source exception branch,
closed launcher and guard PID/start identities, absent scientific outputs and
no matching producer on every system-policy host. The old guard/claim are
preserved with their hashes; GPU locks and foreign processes are untouched.
The case has been resubmitted through the normal queue.

Device discovery also intersects gate availability with one memory inventory
per host before per-slot capacity checks. Unsuitable 16-GB cards therefore do
not delay visits to 32-GB-class V100s/H100s. Atomic GPU reservation, current
RAM promises, disjoint CPU sets, thread caps and scientific gates remain in the
submission path. Exact CPU controllers were stopped only at child-free
boundaries, then resumed; running scientific workers were not stopped.

## Distinct scientific and execution failures

`all_same__occ_4__conf005` was interrupted by a foreign resident during native
work. The original owner is gone and a policy-wide producer audit is clear.
Its exact plan and geometry are continued in a fresh namespace, with the old
partial outputs preserved. A finite independent observer follows that new
attempt; its native G remains pending.

`all_same__occ_2_3__conf006` instead failed a genuine sampled-derivative gate:
1.2138%, above the unchanged 1% limit. CPU component replay rules out the CDS
surface-area and D4 terms as explanations for the approximately 2.88e-6
Hartree/bohr gradient-difference residual. The NumPy mirror of electrostatic
surface retention changes point count at a passing column, not at the failing
column; it is not a verification of GPU electrostatic forces.

Three fresh SCFs at the original base and signed-displacement geometries retain
full solute/solvent/CDS gradient components. Energies reproduce within
7.3e-12 Hartree, but one displaced gradient differs by about 1.995e-6
Hartree/bohr, mainly on the N29--C38 pair (zero-based indices). The repeated
column agrees with the stored Hessian to 8.53e-6 relative error (0.000853%).
This supports a numerical gradient reproducibility problem, while its causal
operator or context dependence remains unresolved. It does not establish an
error in the molecular connectivity or admit the original free energy.

The three-SCF diagnostic completed before an inappropriate Hessian-response
usage assertion failed during context cleanup. The explicit recovery analyzes
only its finished converged SCFs, matching settings and component sums. Its
full execution/GPU-peak qualification remains false. Future gradient-only
diagnostics use the underlying workspace context; the response-use gate for
actual Hessians remains intact. Earlier pre-SCF missing-import and masked-data
failures are also preserved. Preflight now checks the worker's archived module
closure and the complete component-control witness in a read-only scientific
view. No missing control is bypassed.

One fresh full native repeat is prepared from the same checked electronic
parent after this diagnostic, with exact tasks, geometry, numerical settings
and solvent model compared. It must independently pass the normal Hessian,
minimum, displaced-gradient, partition and resource checks. Neither the
diagnostic nor a successful process exit supplies its pending G.

## Sources and reusable commands

Current execution records, immutable selections and raw replay reports are in
`results/validation/polymer_recovery_2026-10-31/`. Scientific calculation data
remain on durable shared storage under:

```text
/nfs/lambda_stor_01/homes/rzhu/redox_polymer_continuations/
  2026-10-31_native_residency_continuation/all_same__occ_4__conf005/
  2026-10-31_bad_column_diagnostic_readonly/
  2026-10-31_native_numerical_repeat/all_same__occ_2_3__conf006/
```

The original native-conveyor failures remain under their original paths. The
three-SCF namespace is a read-only view of its scientific sources, with only
its fresh checkpoints writable. Namespace, plan, owner and source hashes bind
every explicit relocation/recovery.

```bash
python scripts/validation/collect_native_census.py --help
python scripts/validation/diagnose_polymer_components.py --help
python scripts/validation/inspect_pcm_surface_switches.py --help
python -m redox.workflows.polymer_gradient_diagnostic --help
python scripts/validation/analyze_polymer_gradient_diagnostic.py --help
python -m redox.workflows.interrupted_native --help
python scripts/pipeline/prepare_native_numerical_repeat.py --help
python scripts/plotting/polymer_models/plot_occupancy_coverage.py --help
```

M4 remains open. The article's numerical ensemble comparison stays bound to its
existing 75-state snapshot until a complete balanced comparison is ready. The
working census/coverage figure records subsequent verified inputs separately.
