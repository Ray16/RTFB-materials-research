# Maintained workflows and property contracts

Start with `python -m redox workflows` (or `--json`) for a read-only catalogue.
Filter by role with `--category production`, `--category controls`, or
`--category recovery`; `--category historical` isolates compatibility aliases.
Advanced module and script entry points are grouped with their roles and are
listed without importing scientific runtimes or starting jobs.
The catalogue uses the public command registry and fails if a public command has
no declared category. It does not submit calculations.

| Purpose | Entry point | Scientific scope |
|---|---|---|
| Registered molecular campaign | `redox motifs` | Plan by default; execution requires `--execute`, a suitable reserved GPU and state-specific QC |
| Route-admitted primary molecular comparison | `python -m redox.workflows.route_molecular_panel --help` | Exact frozen route/graph/state admission; preserved source geometries and existing corrected electronic worker. A bounded `queue` waits for suitable capacity; own native and reference stages remain mandatory before potential/G publication |
| Route-admitted own native components | `python -m redox.workflows.route_primary_native --help` | Original or independently checked refined electronic parents; one shared constructor, the full 32-GB-class Hessian envelope and unchanged native engine, with actual-center sulfur displacement checks. Failed minima retain blank G; no candidate potential is assembled here |
| Route-admitted geometry refinement | `python -m redox.workflows.route_geometry --help` | Both neutral-singlet saddle branches or the independently diagnosed oxidized-doublet gradient failure; shared strict primary geometry worker/reviewer. Requires 32-GB-class devices and fresh electronic/native work after geometry changes |
| First scoring of accepted refinements | `python -m redox.workflows.route_refined_electronic --help` | Exact neutral-minus/plus and oxidized geometries; same installed/corrected high-level solver, explicit blank historical reference, actual selected occupied-spin review. Geometry energies/nuclear terms cannot transfer; own native checks remain required |
| Resource-controlled shard launch | `redox launch`, `redox run-gpu` | Scheduling/preflight; process completion is not scientific admission |
| Candidate construction | `redox build`, `redox search`, `redox expand` | Graphs and hypotheses; route review precedes new full-property campaigns |
| Exact matched ring panel | `python scripts/candidates/expand_ring_panel.py --help` | Declared baseline/substituents; exact precursor/product and mapped new bond for every graph; retains equivalent recipes, computes no redox properties |
| Unified candidate pool | `python scripts/candidates/merge_registry.py --help` | Fresh source-frozen union with alias deduplication and prospective exact-route decisions; archive and updated route catalog are explicit inputs |
| Table/status review | `redox review`, `redox status`, `redox energy-audit` | Inspect/reconstruct existing outputs without launching QM |
| Low-level calculation stage | `redox uma`, `redox dft`; `python -m redox.qm.sp`, `python -m redox.qm.native_thermal` | Stage controls and defaults differ from the orchestrator; use deliberately |
| Numerical implementation controls | `python -m redox.workflows.screening_control --help` | Explicitly scoped controls; a passed control does not qualify every system |
| Corrected polymer sampling control | `python -m redox.workflows.sampling_energy --help` | Complete 240-input roster; empty/single/double pilots precede same-model energy followups and separate native/ensemble review |
| Corrected polymer native stage | `python -m redox.workflows.sampling_native --help` | Exact checked energy input, full-orbital recovery data, unchanged native physical checks; largest complete native pilot precedes wider native release |
| Corrected native queue | `python scripts/fleet/run_corrected_native_queue.py --help` | Maintained coordinator for the full 240-slot objective; independently collects closed native results, waits for the largest pilot, then prepares and dispatches checked-energy followups |
| Frozen sampling receipt audit | `python -m redox.validation.sampling_snapshot --help` | Verify all 240 declared slots and checked E/native receipts against their source manifest; preserve each native electronic parent and record same-input repeats separately; no fresh physical or ensemble qualification |
| Corrected native batch comparison | `python -m redox.validation.sampling_batches --help` | Analyze complete cells from a source-bound snapshot audit; exact graph/state/thermal context, proper-rotation symmetry RMSD, numerical duplicate and rotor sensitivity, predeclared nested batches. Missing cells and balanced cycles remain unqualified |
| Portable frozen sampling replay | `python examples/occupied_sampling_replay/replay.py` | Uses shipped exact inputs and the shared analysis; repeats graph/RMSD/partition checks without original calculation paths. No new quantum or equilibrium qualification |
| Checked chloride-cohort collection | `python -m redox.validation.contact_cohort --help` | Explicit source-bound starts, signed outcomes and exact native/ledger bindings; group counts come from the checked geometry/energy comparison, never search-hit multiplicity |
| Sampling input export | `python scripts/release/prepare_sampling_inputs.py --help` | Fresh hash-preserving packet from frozen publication inputs; historical producer bytes remain distinct from current reader code |
| Ionic initial-structure replay | `python -m redox.validation.ionic_preparation --help` | Maintained graph, XYZ, stereochemistry and salt-mass audit with an explicit optional original-root relocation. The earlier script remains frozen because historical evidence binds its exact bytes |
| Counterion geometry control | `python -m redox.workflows.counterion_geometry --help` | Oxidized graft plus explicit chloride: +1 singlet, separate from the +1 doublet redox state; geometry/resource pilot precedes wider contacts and native G |
| Counterion small-GPU resource control | `python -m redox.workflows.counterion_gradient_pilot --help` | Prospective exact 42-atom fixed primary SCF/gradient comparison on T4; no optimization, Hessian or native G admission. Full molecular/Hessian defaults retain the 32-GB class |
| Counterion T4 geometry pilot | `python -m redox.workflows.counterion_geometry_resources --help` | One prospective optimization after a checked fixed-gradient witness; other prepared repairs require closed full geometry resource/physical review. No Hessian or native G admission |
| Counterion thermochemistry | `python -m redox.workflows.counterion_thermochemistry --help` | Prepare/launch/queue own matched primary electronic/native components after the full pilot; shared owner-aware bounded dispatch and unchanged 32-GB class. Analyze only complete independent reviews, with explicit two-to-one solute conversion |
| Saved ionic Hessian recovery | `python -m redox.workflows.charged_hessian_recovery --help`; `python scripts/fleet/run_saved_ionic_validation.py --help` | Recovery utility for one proven closed post-Hessian file-lock failure; retain exact full response/Hessian, require a fresh center and every derivative/minimum check. Uses the common probe/launcher and the full declared ionic resource class |
| Signed minimum repair | `python -m redox.workflows.counterion_mode_repair --help` | Independently diagnose a closed physical failure and prepare both signed-mode branches; fresh geometry and native QC remain mandatory |
| Ionic motif construction | `python -m redox.workflows.viologen_motif_fragments --help` | Reuse the audited methyl-cap constructor at four checked parent frames; exact CH2 attachment, retained coordinates and charges are checked |
| Frozen motif electronic control | `python -m redox.workflows.motif_interaction_energy --help` | Prepare/launch/collect matched A/B/AB gas and MeCN electronic triplets after a checked pilot; no binding G or complete backbone decomposition |
| Charged-polymer pilot | `python -m redox.workflows.charged_polymer_pilot --help` | Route-bound private namespace, first-class resource and physical checks |
| Directional solvent response | `python -m redox.workflows.solution_pair_namespace --help` | Exact matched candidate/state pair, independently reviewed |
| Historical compatibility entry | `python -m redox.current` | Alias of the same public CLI, not a second production implementation |

Historical source snapshots and recovery namespaces remain intact. Inspect exact
host/PID/start identity, checkpoints and source contracts before recovery. A
prepared namespace, launch receipt or old status file does not establish a live
calculation, an admitted G or an ensemble result. Only expand a resource class
after its declared pilot passes. The unified design pool is
`results/candidates/pool/current.json`; a route pass is distinct from a full
material decision.

The route-admitted molecular queue and saved-ionic validation coordinator share
`prepared_namespace_queue`: exact owner closure for resume, bounded resource
refusal accounting, source guards, shared storage/GPU admission and preservation
of every started scientific attempt. A multi-input queue exposes only its next
pending namespace to the advisory priority rule, so it cannot hold extra slots
for already running states. Scientific admission remains in each caller.
Equivalent registry declarations may differ in SMILES spelling or source label;
their exact structures and charge/spin/state hypotheses must all agree before a
campaign starts. Conflicting declarations require review.

Independent states continue to observe or dispatch when a partner has failed;
the failed attempt is never restarted by this rule. A terminal incomplete queue
is reported once all remaining states have finished or require explicit review.
Added finite-difference probes must use the native engine's recorded center
coordinates, preserving its pinned unit convention rather than applying a
second Angstrom-to-Bohr constant.

The shared remote launcher has an explicit adapter for the maintained route
geometry entry point; existing rows retain their original screening-control
launcher by default. Each adapter enforces its own chemistry/state/resource
contract before calling the same reservation and budget machinery. It cannot
select an arbitrary Python module.

## Physical context at arithmetic boundaries

`redox.properties.context.combine_properties` is the shared checked-arithmetic
boundary introduced for new property combinations. It retains the component
identity, charge/multiplicity, source bindings and QC together with the quantity,
unit, electronic protocol, thermal definition, solvent/model, temperature,
standard state and reference. Mixed contexts, nonfinite values, missing state
information and pending QC are rejected. Explicit conversions must precede
combination; the helper never invents a conversion or validates a raw calculation.

The checked quaternization publisher now uses this contract at every cutoff and
for every precursor seed. It first verifies the actual quantum records and only
then declares independently checked component inputs. Its result retains those
contexts and references the verified source inventory. The legacy pure
`quaternization_delta_g` arithmetic API remains usable without context records,
but explicitly reports `property_context_verified=False`; it cannot acquire
physical qualification merely through arithmetic. This is an initial integration,
not a claim that all historical property readers have migrated to the contract.
Gas diagnostics, directional solution relaxation, composite free energies,
conditional occupancy differences and electronic couplings remain different
quantities. Never silently add a solvent contribution to an already composite
solution-relaxation value or promote an electronic occupancy difference into G.

## Runtime discovery and portability

The common resource layer discovers the system reservation executable from an
explicit `REDOX_GPU_RESERVE`, then installed site defaults in `config/runtime_paths.json`, then `PATH` and the account's `~/bin/gpu_reserve`.
An explicitly configured invalid tool fails rather than falling back. CPU-only
imports work when no GPU service is installed; GPU execution still cannot bypass
the missing gate. `GPU_LOCK_ROOT` is the shared ownership-namespace override
consumed by both the Python owner check and the system gate. The default follows the matching site lock directory when that site gate is selected, otherwise `~/.gpu_locks`. All workers must inherit the same authoritative tool/config and
lock namespace. This is runtime discovery, not permission to replace the system
reservation policy or alter another job's ownership.

The generic resource and full-molecular launch paths use these declarations.
Frozen scientific namespaces can retain the original cluster paths for evidence
integrity; their explicit relocation and pinned-runtime contracts are separate.
Dependency/environment qualification and architecture-specific numerical tests
remain necessary on a new host. Numerical and publication builds must not rely on
GPU availability when their task is CPU-only.

The charged-native followup entry point requires the exact independently checked
pilot, complete G values, matching physical context and identical worker/kernel
bytes, and admits only dimensions within that pilot. It reuses the established
node preflight/launcher, capacity utility and storage accounting. A RAM refusal
before worker start is not a scientific attempt. The unstarted-namespace
relocation entry point preserves the original and accepts only hash-identical
dependency replacements; an existing scientific start prohibits this operation.
This permits shared storage to replace a node-local package location without
changing the method, geometry, state, model, or original evidence.


The exact 137-atom, charge-4 double-graft native class uses an explicit separate
H100 dimension pilot (`charged_polymer_pilot --native-from-review ...
--larger-dimension-pilot-review ... --out ...`). Its state/dimension, retained
output, device and peak limits are shared by preparation, launcher, frozen
worker helper and independent review. The completed smaller native witness is
required to prepare it; the *larger pilot itself* must pass before matching
larger followups may dispatch. The same followup entry point handles each
class through a source-bound manifest, with no automatic smaller-to-larger
admission. The existing compact numerical adapter and physical QC are unchanged.
Review inputs are copied into each new namespace to stay visible inside the
isolated source checkout.


## Matched native analysis and primary implementation controls

`python -m redox.workflows.ionic_native_cycles --help` is a read-only analysis
entry point. It verifies exact graphs, charges, model, implementation, thermal
protocol and source receipts before context-checked sums. Missing partners keep
the whole cycle blank. Finite-minimum G is distinct from an occupied ensemble,
isolated motif binding or electronic coupling.

`python -m redox.workflows.screening_native_panel prepare --help` prepares
source-bound primary molecular repeats; `dispatch` requires the exact largest
closed/open-shell physical/resource pilot reviews. Both this panel and charged
native followups reuse one bounded launcher, storage accounting, authoritative
reservation queries, node preflight and CPU/RAM guards. A slower healthy-node
query may be retried once through the same system gate; an observation timeout
never reclaims a worker or changes a reservation.

Actual SCF molecular frames, occupied electron counts/global S², execution
closure and resource identity use shared CPU reviewers. Their as-used source
bytes are frozen beside each receipt. The independent full-precision spin gate
is separate from comparison with four-decimal producer diagnostics.

The large RKS recovery uses `charged_partial_recovery` after a proven closed
H1 workspace failure. It restores exact orbitals/grids and the completed solute
partial tensor, while leaving the full solvent H1 method intact. A prospective
16-GiB RKS workspace and 40-GiB compact budget require the declared H100 class;
all derivative, residual, symmetry and minimum limits remain unchanged.

Quantum preparation discovers a declared shared GPU4PySCF source root from
`REDOX_GPU4PYSCF_ROOT` or site defaults, falling back to the actual installed
package only when none is declared. Grid, integrator and compiled-kernel bytes
must match before shared-path binding. The discovery source/config is frozen
in the namespace; this does not change the executable or numerical method.

## Corrected sampling handoff

The electronic and native coordinators share the append-only storage-peer
registry. The native coordinator consumes the complete electronic ledger while
retaining every pending slot. A checked electronic input must match the original
task, geometry, model and protocol, including when its execution namespace was
recovered. A completed largest native pilot and identical calculation helper
bytes are required before preparing the remaining native inputs.

`redox.validation.sampling_native_receipts.checked_native` is the arithmetic
reader for this stage. It checks both electronic and own nuclear components,
all three rotor settings, exact geometry, charge/spin, solvent/model, electronic
and thermal definitions, temperature, translation convention and field-specific
QC. The collector retains that context alongside each G. Pending G stays null;
240 native completions still require separate duplicate, sampling-convergence
and balanced ensemble analysis.

The queue reuses the shared reservation/preflight dispatcher and its RAM/storage
accounting. Unknown/live owners and scientific failures are preserved. Only
proven prelaunch refusals can retry, with cooldown and a three-receipt limit;
three terminal failures block the affected host. Source changes stop further
admission while preserving workers. Explicit resume requires the exact previous
controller to be proven stopped. All time windows and paths are declared CLI
arguments, and a new controller revalidates receipt hashes before caching them.
Already prepared native inputs dispatch before additional preparation. Source
checks run again before preparation and dispatch; the status records which
phase is active. A slow preparation pass does not change the scientific
qualification of already running calculations.

The counterion geometry entry also supplies `collect`, which retains all seven
starts and field-specific pending values. Geometry-level E differences are not
association G. Thermochemistry reuses the existing molecular SCF and full native
Hessian/partition routines, with shared actual-checkpoint and bounded-launch
adapters. The first full contact component gates further thermochemistry.
`analyze` checks atom/charge balance, exact graft retention, methods, thermal
definitions, numerical implementation and physical-review backing. A bound
cluster is one statistical solute despite its disconnected SMILES components;
the 1-atm-to-1-M correction is applied after the balanced source-convention sum.
Each starting pose remains a separate finite scenario, without inferred
experimental populations, graft conversion or yield.

The counterion failed-mode workflow preserves the original SCFs and Hessian,
reconstructs their sampled derivatives and signed internal modes, and retains
both ±0.1-Å branches. The repair uses stricter optimization convergence without
changing the molecular Hamiltonian or the −5-cm⁻¹ native minimum gate. Stationary
geometry, native minimum and composite G retain separate qualification fields.

The full counterion pilot measured approximately 1.65 GiB host RSS and 20.9 GiB
GPU usage. Matching followups may use a declared 16-GiB host envelope after
checking 25% headroom against the measured pilot, while retaining the validated
32-GB-or-larger GPU requirement. These resource bounds do not alter the SCF,
Hessian, finite-difference or thermal protocols.

## Prioritizing a pending lead

`config/campaign_dispatch_priorities.json` declares an exact prepared namespace
and its controller status. The shared dispatcher may defer one compatible free
slot to a verified live pending lead. It checks memory/architecture, CPU/RAM and
aggregate storage, and continues allocating other suitable slots. A started or
inactive priority cannot hold another GPU. Each pass freezes the configuration
and the as-used capacity, ownership, launch and storage helpers. No site GPU
policy, reservation ownership or scientific method is changed. The shared prelaunch accountant separates explicit closed GPU/RAM/storage
refusals from bounded environment failures for electronic, native and saved-ionic
controllers. It checks frozen source branches before Popen, exact owner closure
and absent scientific outputs. Resource waits preserve cooldowns without
consuming scientific/prelaunch-failure budgets; uncertain starts require review.
