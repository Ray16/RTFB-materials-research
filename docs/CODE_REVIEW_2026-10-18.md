# Code review and refactoring

## Current checkpoint: 2026-10-19

The full software suite passes **514 tests**. The latest additional fix makes
alternate-spin helper repairs retain both optimizer and final-SCF convergence,
refuse replacement geometries on either failure, and clear inherited direct-energy
producer declarations. Five regression cases cover successful convergence and
both failure paths, including preservation of existing helper files. The previous
producer is retained in
[the source-epoch record](../results/validation/spin_helper_convergence_2026-10-19/source_epoch.json).
Helper repair adoption also reuses the qualified transaction installer, carrying
`projected_thermal` with the SP/native caches and preserving old bytes on failed
validation or changed input. Five further regression cases cover that adapter.
All primary parent workers passed their initial source guards before these changes;
their geometry/energy/thermal numerical producers were unchanged.

The maintained pipeline now binds energies, native partitions and equilibrium spin
policies to actual inputs, validates complete replacements before atomic promotion,
and preserves optimizer convergence and original records. Both density-cache
schemas require explicit provenance. Shared publication inventories, CSV/workbook
schemas and readiness checks preserve missing fields and field-specific QC.

The strict experimental loader now supports positive-electrode oxidation ladders
and binds experimental CSVs into table-review fingerprints. The new primary
phenothiazine source supplies five strict target potentials from four parents;
irreversible waves and conflicting values remain excluded from strict scoring.
The [parent panel](../results/validation/phenothiazine_parent_panel_2026-10-19/README.md)
records its calculations and independent partial reviews. Started calculations do
not establish matching accuracy or the preferred spin surface.

The last full table snapshot passes **30,455 arithmetic checks**, and the additional
candidate export passes **572 checks**, with zero discrepancies. These are scoped
snapshots; new parent calculations require a fresh assessment after their QC.
[Current evidence](../results/LATEST_REVIEW.json) identifies the exact receipts,
24/25 admitted additional potential pairs, source catalog and checked slide deck.
Both remaining hexyl-NDI branches failed physical gates. Host-only diagnostics
reproduce the problem in the SMD CDS term; neither geometry is admitted.
Solution kinetics, broader experimental accuracy and exact-graft validation remain
unfinished. User-staged archive changes remain separate and untouched.

## Gas-reference solver correction and remaining forwarding

Three gas-reference workers now use a shared strict solver contract: finite SCF
energies, optimizer and final-SCF convergence, and full grid-response gradients.
They preserve incompatible/legacy caches, bind reuse to method and source identity,
store all four energies and geometries, and fail incomplete CLI jobs. The generic
GPU HOMO reader now transfers device arrays explicitly and withholds usable
reorganization values for unknown or unbound anions. The acene protocol records
its actual requested basis. An actual CPU solver control and 48 new regression
cases pass; see the [qualification](../results/validation/gas_reference_contract_2026-10-19/README.md).

The read-only inventory covers 482 saved reference records without complete new
solver verification. This is insufficient evidence to certify them, not proof that
their numerical values are wrong. Existing datasets are unchanged. The public
`properties/nelsen.py` forwarding remains pending its active source epoch; its
production consumer uses only `embed`. The source inventory and resolution are
kept under [remaining solver review](../results/reviews/2026-10-19/remaining_solver_review/).

## Earlier review checkpoints

The dated material below retains the original findings, investigations and counts.


## Current corrections: 2026-10-19

CN-PMDI now has a selected, source-bound two-spin equilibrium policy, verified by
both real readers. The triplet population is 32.16% at 298.15 K; the mixture uses
each component's native spin degeneracy once. Every included component pair passes
solution binding and distance-graph integrity. Kinetic values remain flagged.
Both ongoing methyl-NDI branches passed their source guards before the shared
analysis hook changed. A follow-up export test caught first-row-only CSV/workbook
schemas; both now retain all optional QC fields. The full suite passes **429 tests**.
See the [equilibrium implementation and proof](../results/validation/spin_ensemble_integration_2026-10-19/README.md).
Earlier paragraphs below preserve the preceding checkpoints.


Near-spin native surfaces previously stayed blocked because primary free energies
had no equilibrium selection. A source-bound two-spin policy now supplies the
mixture correction to both readers and ranking gates, with transactional selection
and stale-input refusal. Component spin degeneracy is counted once. The current
CN-PMDI policy was independently replayed and selected; kinetics remain withheld.
Integrity now checks all included spin-component pairs. The independent audit
also reconstructs distance bond graphs instead of trusting table bond counts.
The full suite passes **428 tests**. Earlier checkpoints below retain their counts.


The generic native repair path also omitted `projected_thermal` when promoting a
replacement and changed canonical geometry before final admission. It now prepares
and validates the complete replacement, retains old addressed histories, swaps
state directories with rollback, and preserves every archived byte. Repair
optimization uses the kernel convergence result and reoptimizes resumed checkpoints
rather than declaring them converged. Four promotion regression cases cover the
complete partition, pre-install refusal and post-install rollback. The full suite
now passes **395 tests**. Two private NDI radical minima independently pass 157
checks each and have been promoted with original-state archives.


Native thermochemistry previously included rigid-motion modes. The corrected
internal-mode partition is integrated into the producer and both primary readers,
with 232 cases passing 6,961 independent checks. The installed xTB thermodynamic
kernel runs in an isolated CPU subprocess, avoiding the PySCF/D4 shared-library
conflict while retaining genuine soft internal modes. Raw evidence is preserved.
See the [qualified thermal correction](../results/validation/rigid_mode_thermal_2026-10-18/README.md).

The geometry driver also passed a mean-field object to geomeTRIC with the installed
gradient default `grid_response=False`, then used `optimize()` which discarded
optimizer convergence. It now passes an explicit full-response gradient scanner to
`kernel()`, records geometry and SCF convergence, and refuses unconverged strict
results before overwriting geometry. Seven maintained regression cases and an
actual CPU water optimization cover this correction. The full suite passes
**391 tests**. Published native quantities retain independent stationary-gradient,
minimum, derivative and spin gates; historical geometry records are not relabeled.

The candidate feed now admits 22 of 25 complete potential pairs under corrected
thermochemistry. Remaining native minima, close-spin decisions, solution kinetics,
and exact-graft accuracy are unfinished. The corrected global arithmetic review
passes 28,967 checks with stable inputs and zero discrepancies. Later
sections preserve earlier review checkpoints and their historical counts.

The subsequent [scientific continuation](SCIENTIFIC_PROGRESS_2026-10-18.md)
records the implemented benchmark correction, live SP revalidation and six new
matched-comparison candidates. The review below preserves the initial findings
and first refactoring checkpoint.

A later source-lineage review found 92 generic repair records pointing back to
their replaced canonical geometry. All 92 originals were recoverable from
explicit, byte-verified archives. The promotion helper now writes the archive
path, canonical `opt_out`, helper-source snapshot and atomic metadata; it does not
copy old direct-energy producer declarations into repaired records. Historical
bytes remain unchanged, and `redox status --include-lineage` resolves their
original geometries by hash. The suite passes 332 tests. Numeric native/SP
producer modules and physical thresholds were unchanged by this metadata fix.

A subsequent scorecard review found that unfinished or stale alternate-spin helpers
could override UMA confidence using their electronic energies. Candidate exports
already rejected these helpers. The scorecard now requires a completed current
source-bound comparison, matching helper protocol/multiplicity, energy and native
thermal dependencies, with the same near-gap publication gates. Metadata predicates
are shared in `validation/spin_evidence.py`; the raw audit reconstructs identities
and gaps separately. Live kernel writer locks also block ranking while helper
directories are being replaced. All 140 registered molecules agree between readers; four
unresolved legacy comparisons are excluded from primary ranking. The full suite
passes **343 tests**. Progress labels now distinguish partial native, active-SP and
geometry records rather than reporting every partial result as gas UMA only.

Reviewed the working tree on 2026-10-18, including existing uncommitted work.
The baseline contains 41 package modules, 101 Python scripts, 17 test files,
and 10 Python configuration files. All 181 baseline tests pass in the existing
NFS `redox` environment with native thread pools limited to one thread.

The main maintenance problem is that operational workflows have grown in
scripts while shared data ownership and supported entry points remain implicit.
The existing protocol-addressed energies, explicit incomplete values, separate
spin and thermal gates, and independent arithmetic audit are valuable contracts
to preserve during cleanup.

## Findings before edits

| Priority | Finding | Trigger and consequence | Correction |
| --- | --- | --- | --- |
| High | Single-point adoption lacks source identity checks | `qm/sp.py:ensure` checks method metadata but does not bind the adopted `result.json` energy to the current geometry or requested charge/multiplicity. Replacing `opt.xyz`, or requesting another spin, can label an old energy with a new address. | Require producer-recorded geometry and full energy protocol, matching charge/spin, finite energies and explicit convergence. Preserve existing addressed records. |
| High | Shared manifest lock is used by only one builder | `build/candidates.py` locks `.manifest.lock`; `build/standalone.py` and `build/validation.py` perform unlocked read/merge/truncate writes. Concurrent builders can lose entries or expose a partial manifest. | Share the same lock and atomic CSV replacement across all three builders, holding the lock through read/build/merge. |
| High | Finalization can start an unreserved GPU calculation | `scripts/pipeline/finalize_after_dft.sh` invokes `reorg --all-cross --backend gpu` directly. The default command can start compute without the reservation and node preflight required by `AGENTS.md`. | Make finalization consume existing cross points; route new calculations through the supported reservation-aware workflow. |
| Medium | Isolated table reviews leave mutated module globals | `scripts/pipeline/review_tables.py` changes imported modules' `RESULTS` paths permanently. Later calls in the same interpreter can write to an earlier review, and concurrent calls race. | Put table orchestration in the package and run stages in separate processes with explicit destinations. |
| Medium | Record and run-manifest writes can be observed halfway through | SP, cross-point and run-manifest writers use `write_text` on live JSON files. A watcher or interrupted process can see truncated evidence. Several other scripts implement their own temporary-file pattern. | Introduce one atomic I/O utility and adopt it on these production paths. Atomic replacement alone does not establish exclusive writer ownership. |
| Medium | Alternate SP protocols can claim solver settings that were not applied | `sp.ensure(proto=...)` forwards functional/basis/solvent, while density fitting, SCF tolerance, solvent model and iteration limit remain hardcoded in `dft_smd`. | Reject unsupported solver overrides before writing or computing. Record the actual producer protocol. |
| Medium | Older launchers do not implement the complete current resource contract | Thread limits in `run_dft.sh`, `run_uma.sh`, and fleet shell launchers do not establish disjoint CPU affinity, RAM suitability or validated GPU memory. Reservation alone does not prove node suitability. | Identify supported launch paths and retire or harden legacy shortcuts before future use. Do not validate a launcher by running an unrequested campaign. |
| Medium | Documentation disagrees with the working tree | README/map claim bulk `calcs/` is ignored and `archive/` preserves old code; Git tracks hundreds of calculation artifacts and the current archive is empty. The map says scripts contain no reusable logic, yet scripts contain about twice the package's Python lines. | Describe actual ownership and navigation, mark historical documents, and move reusable orchestration into the package incrementally. Preserve calculation paths and evidence hashes. |
| Medium | Environment setup is not reproducible as described | README/setup comments advertise CUDA 12.4 while setup defaults to `cu128`; many runtime libraries are unpinned and the validated thermal protocol requires specific PySCF/GPU4PySCF versions. | Align setup instructions and separate a recorded, validated environment from best-effort installation. Preserve the current environment while correcting setup. |

## Review coverage and limits

Inspected structure construction, shared readers, active SP protocol and adoption,
DFT batch selection, thermochemistry admission, gas and solution reorganization
interfaces, potential/disproportionation aggregation, scorecard/Pareto contracts,
verified candidate collection, publication freshness, table review, shell launchers,
environment setup and repository guidance. Parsed all Python files. The tests cover
thermodynamic identities, malformed/stale caches, molecular identity, spin admission,
derivative evidence, complete ranking objectives and publication withdrawal.

Passing software checks does not establish experimental accuracy or a physical
minimum for every calculation. Correctness work must distinguish arithmetic,
source/protocol consistency, physical QC, and validation against measurements.
Keep missing quantities blank and retain parent measurements as separate evidence.

## Refactoring boundaries

`src/redox/` owns reusable implementations; `scripts/` provides operational entry
points and historical campaign drivers. Configuration declares chemistry and
protocols. `library/` and `calcs/` are scientific evidence with existing consumers;
`results/` contains derived tables, immutable snapshots and campaign evidence.
Do not relocate those trees or broadly alter Git tracking during a code refactor.
The independent audit must retain its own arithmetic and thermal admission checks.

The initial refactor addresses shared I/O and manifests, SP identity, isolated table
orchestration and navigation. Subsequent work audits present results, fixes remaining
launch/setup issues, and defines candidate expansion using exact graft identities,
coverage and resource requirements. Completion and validation are recorded below.

## Validation after refactoring

The refactor adds a discoverable `python -m redox --help`/`redox` CLI, package
workflows for isolated table review and resource-controlled execution, shared
manifest transactions and atomic record I/O. Geometry mapping now lives in
`build/geometry.py` with compatibility imports in `core/common.py`.

All three builders use the same manifest lock. New SP adoption requires a
producer-recorded geometry hash, complete energy protocol, matching charge/spin,
finite energies and explicit convergence. Unsupported solver settings are
rejected, and a geometry changing during an SP or cross-point calculation cannot
publish an energy under its earlier address. Existing valid addressed records
remain readable.

GPU discovery delegates to the system gate. DFT/UMA and single-command launchers
share CPU allocation, RAM/GPU suitability checks, thread limits and ownership
records. Allocation/start is serialized locally; the system gate remains the
atomic GPU owner. GPU library probes execute inside the reservation. Historical
claim-based shell queues and the unreserved OROP fanout are retired; their source
is preserved as reference text in `scripts/legacy/`. Full molecular jobs require
the proven 32-GB-class devices. Candidate preparation no longer imposes its old
one/two-worker cap. The historical single-node deferred controller remains a
dated campaign implementation, not the default expansion policy.

The environment bootstrap now agrees with its actual CUDA-12.8 default and pins
the native implementation's PySCF/GPU4PySCF versions. Development and navigation
guides distinguish maintained code, compatibility wrappers, scientific evidence
and historical campaign documents. Existing calculation paths, records, staged
changes and manuscript/slides source are preserved.

The final CPU suite passes 210 tests. A real adapter smoke test on lambda2's V100-32GB passed dependency checks,
ran with exactly CPUs 0–1 and two native threads under `gpu_reserve`, and exited
successfully. It performed an environment probe rather than a quantum calculation.
See `results/reviews/2026-10-18/gpu_adapter_smoke/` for execution evidence.

## Current property audit

The new table workflow rebuilt every derived table and workbook into
`results/reviews/2026-10-18/tables_corrected/`. The independent audit passed all
**19,753 checks with zero discrepancies**. The first review found two display
differences at the PMDI mass value 402.45: the audit applied one-decimal rounding
to an unrounded binary sum, whereas the scorecard consumes a four-decimal CSV.
The independent reconstruction now follows that declared serialization precision;
the underlying capacity formula and scientific masses were unchanged.

`results/reviews/2026-10-18/readiness.csv` reports 315 declared states: **292 ready,
21 incomplete, two invalid**. The incomplete states lack admitted native thermal
records. The invalid controls are `mv_ip2/s0` (wrong selected spin) and `mv_ip1/s0`
(changed connectivity); these ion-pair helpers are excluded from candidate
capacity/scorecard screening. The audit also lists 113 record/physical annotations:
73 gas-unbound anions, 12 neutral-solvation outliers, 21 missing thermal terms,
four missing active records in the wider on-disk tree, two monoion-solvation
outliers, and the reported Fc gas IP. These have different meanings and do not
justify replacing missing values or changing computed energies.

Arithmetic agreement is established for this snapshot. Physical admission,
experimental accuracy and complete candidate coverage remain distinct. A further
scientific pass must also review the historical OROP benchmark's use of legacy
GFN2 corrections against a stored Fc reference now maintained under native
thermochemistry; the old benchmark file has no thermal-definition metadata.
There are also 143 historical adopted SP files in the broader DFT tree. New
adoptions now require source geometry proof; older addressed records retain their
existing provenance and need source-level review or recomputation where that proof
cannot be recovered. That benchmark file is a historical input in this isolated
review, not new native benchmark validation. Recheck benchmark/reference matching before regenerating it or using
it as current calibration evidence.
