# Polymer execution profiles and checked state handoffs

The standard polymer optimizer, electronic-state reader, structural audit and
native thermochemistry retain their existing acceptance criteria. Execution
profiles control SMD allocation order, workspaces and bounded contractions.
Small implementation controls and an exact-graph resource point establish
different evidence from a stationary state or a converged conformer ensemble.

Use [workspace_polymer_campaign.py](../scripts/pipeline/workspace_polymer_campaign.py)
to continue a completed bounded-gradient resource point in fresh calculation
storage:

```bash
python scripts/pipeline/workspace_polymer_campaign.py prepare-profile \
  --source CHECKED_RESOURCE_CAMPAIGN --out NEW_FULL_CAMPAIGN
```

The source point must be finished, within its declared memory budget, bound to
the checked adapter and backed by both spin controls. The new campaign retains
the original task, geometry, physical model and optimizer. Its copied pilot
receipt explicitly identifies the original process and states that no new
pilot was performed. Full geometry-force storage is separate, so the existing
pilot evaluations cannot collide with optimization output.

A separately completed electronic workflow can supply a native thermal stage:

```bash
python scripts/pipeline/workspace_polymer_campaign.py prepare-native \
  --source CHECKED_ENERGY_CAMPAIGN --out NEW_NATIVE_CAMPAIGN \
  --controls CHECKED_EXECUTION_CONTROLS --workspace-GiB 64 \
  --minimum-GPU-MiB 81920
```

This checks the stationary electronic report, structural/force evidence, active
energy on the exact saved geometry, backing hashes and completed workflow status.
The new energy receipt points to the original checked report. Optimization is
not repeated and no G is created by the handoff. The normal native resource
pilot, DFT Hessian, finite-gradient checks, projected modes and partition gates
remain necessary. The current adapter is scoped to one neutral-singlet graph
per campaign; it is not a qualification of charged polymer models.

Both prepared campaign types use the existing `launch`, `worker` and `collect`
commands and shared GPU reservation mechanism. Host admission uses the measured
density-fitting storage estimate for electronic work and the declared 256 GiB
native-thermal budget. Actual device, CPU affinity, native thread limits and
process ownership remain part of launch verification.

The [candidate observer](../scripts/pipeline/coordinate_candidate_baseline.py)
can declare an external energy continuation with a native target, checked
controls and execution budget in `baseline_recoveries.json`. It registers the
thermal campaign only after the source workflow completes its checks. Existing
targets must retain the same task and numerical model; partial or mismatched
targets are preserved for review. Active or unresolved task claims suppress
extra workers. Source evidence and historical registry versions remain intact.

Current concrete branches are indexed from [results/README.md](../results/README.md).
Passed point/resource controls, checked E, checked native G, conformer convergence
and synthesis admission are separate statuses. Missing quantities remain blank.

## XC temporary lifetimes and measured-cavity native admission

The [XC campaign wrapper](../scripts/pipeline/xc_workspace_campaign.py) composes
the checked solvent-workspace profile with a process-local release of the final
AO batch after the completed density loop. It records transformed execution
bodies and their hashes; no installed library is changed. Separate original/
modified spin controls and sampled checkpoint total-energy derivatives are
required before using a saved failed trajectory as a warm start. The preparation
entry point is `scripts/pipeline/prepare_xc_checkpoint_resume.py`. Only starting
coordinates are reused: optimization and all ordinary physical checks repeat.

For native thermochemistry, this wrapper can replace blanket 256 GiB admission
with the original bounded-PCM forecast evaluated at the parent's measured maximum
cavity size plus 1,024 safety points. Cache and headroom limits remain 128 GiB
each. Both effective and legacy estimates are recorded. The site-3 case forecasts
216.16 GiB; hardware suitability still requires an actual native resource pilot.
These controls do not qualify charged polymers or every possible cavity.

The candidate registry supports per-case host lists and checked-energy handoff
from an external continuation. It only adopts a native target after exact-task,
model and completed-parent checks pass. Existing live/unresolved claims suppress
duplicate workers. Source-pin checks remain active for both drivers and the
observer. Running producers must not be edited to retrofit a recovery.

## Explained coarse-step force failures

`qualify_saved_polymer_energy.py` creates a separate, source-bound energy capsule only when the coarse discrepancy is quantitatively attributed to the same CDS energy term and a predeclared independent full-quantum refinement passes. Both original components require four accepted displacements (0.0005, 0.0001, 0.00005 and 0.00001 bohr), with unchanged 1e-6 Ha/bohr error tolerance. Baseline repeatability and finest-step consistency add checks. The original failed 0.001-bohr result is retained explicitly. Changed geometry, unexplained residuals, missing coverage or failed fine checks cannot be rescued by this procedure.

The capsule changes the force-validation procedure, not the saved geometry or electronic Hamiltonian. Its receipt names this policy; it makes no Hessian-minimum or G claim. Native thermochemistry retains the existing separate physical gates and output storage. Historical primary evidence is unchanged.

Preparation subprocesses recover the observer's initial permitted CPU mask before
their normal allocation guard selects a separate free CPU. This avoids trapping
the child on its pinned parent's occupied CPU; scientific workers keep their
own fixed allocations. Resuming a candidate observer retains its recorded
finite deadline exactly, rather than extending the observation window.

For an imported energy receipt, run
`materialize_imported_energy_parent.py NATIVE_CAMPAIGN` after preparation to make
the native report's canonical parent path resolvable. It copies the already
checked report byte-for-byte, refuses conflicting destinations, retains the
original receipt/source and records both paths and hashes. It changes no plan,
scientific result or running worker. Both current single-site native campaigns
have these checked aliases. Integrating this step into the shared preparation
driver remains deferred while that producer is source-pinned by live native jobs.

## Bounded first response and exclusive continuation

The new `first_response_polymer_campaign.py` profile requires both full spin controls and an independently reconstructed fixed-density response test for the exact saved geometry. It composes bounded rectangles from the unchanged installed pairwise primitive with contracted Hessian R applications. It materializes imported parent aliases automatically. Complete native resource, Hessian, minimum and partition gates remain separate; the response component is not a G certificate. The isolated response test uses its declared physical device ceiling, while the native pilot retains its own stricter resource criterion.

`watch_gpu_exclusivity.py` adds finite CPU-only post-launch residency audits to the shared allocation gate. Delayed foreign residency stops only the owned matching PID/start/campaign, leaving peers, claims and reservations untouched. Candidate phases with interrupted workers remain held until an explicit fresh continuation is prepared. `prepare_exclusive_checkpoint_resume.py` verifies termination, preserves the old trajectory, rechecks checkpoint connectivity/stereochemistry and binds a live residency observer. Its warm-start coordinates remain unqualified; previous same-graph resource evidence is not relabeled as an exact-checkpoint energy calculation. All final gates repeat normally.

## Later switch curvature and full-native qualification

`switch_curvature_polymer_campaign.py` additionally uses unchanged all-atom primitive evaluations in grid tiles and host-backed switch tensors. Native curvature bindings are refreshed at execution time so the later certified bounded-native context is honored. Full guarded spin/Hessian controls must exercise host storage, live bindings and a durable solute-component checkpoint. The declared native measured-peak ceiling is separate from the workspace estimate and preserves at least one GiB above it in the minimum physical device class. This resource limit is not a relaxation of native physical criteria.

A passed response-only test cannot qualify later arrays: the exact site-3 geometry requires 43.37 GiB for its two dense switch tensors. Their host-backed successor still needs a complete large native pilot. Partial solute checkpoints remain unqualified until full assembly and normal checks finish.
