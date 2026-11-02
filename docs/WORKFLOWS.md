# Maintained workflows and property contracts

Start with `python -m redox workflows` (or `--json`) for a read-only catalogue.
The catalogue uses the public command registry and fails if a public command has
no declared category. It does not submit calculations.

| Purpose | Entry point | Scientific scope |
|---|---|---|
| Registered molecular campaign | `redox motifs` | Plan by default; execution requires `--execute`, a suitable reserved GPU and state-specific QC |
| Resource-controlled shard launch | `redox launch`, `redox run-gpu` | Scheduling/preflight; process completion is not scientific admission |
| Candidate construction | `redox build`, `redox search`, `redox expand` | Graphs and hypotheses; route review precedes new full-property campaigns |
| Table/status review | `redox review`, `redox status`, `redox energy-audit` | Inspect/reconstruct existing outputs without launching QM |
| Low-level calculation stage | `redox uma`, `redox dft`; `python -m redox.qm.sp`, `python -m redox.qm.native_thermal` | Stage controls and defaults differ from the orchestrator; use deliberately |
| Numerical implementation controls | `python -m redox.workflows.screening_control --help` | Explicitly scoped controls; a passed control does not qualify every system |
| Corrected polymer sampling control | `python -m redox.workflows.sampling_energy --help` | Complete 240-input roster; empty/single/double pilots precede same-model energy followups and separate native/ensemble review |
| Corrected polymer native stage | `python -m redox.workflows.sampling_native --help` | Exact checked energy input, full-orbital recovery data, unchanged native physical checks; largest complete native pilot precedes wider native release |
| Counterion geometry control | `python -m redox.workflows.counterion_geometry --help` | Oxidized graft plus explicit chloride: +1 singlet, separate from the +1 doublet redox state; geometry/resource pilot precedes wider contacts and native G |
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
