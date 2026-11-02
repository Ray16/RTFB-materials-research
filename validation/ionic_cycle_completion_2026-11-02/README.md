# Exact lead thermochemistry and occupancy completion

The [current eleven-component collection](partial_cycles_11.json) closes three
balanced finite ionic viologen G cycles. At the 50-cm⁻¹ rotor setting they are
−3.9365 kJ/mol (all-same adjacent), +22.5790 (all-same one spacer), and −5.2972
(alternating one spacer). The 25/50/100 values are preserved in the collection.
The all-same adjacent/spacer electronic contributions are −3.4607/+18.6848 kJ/mol.
These whole-chain responses include deformation and solvent changes; they are
not isolated motif binding, coupling, converged ensembles or material admission.

All three 137-atom +4/singlet double-graft components pass full Hessian, minimum,
sampled derivative, own nuclear-term and all 25 actual SCF-frame checks.
Their minimum frequencies are 6.9604, 6.6572 and 8.4292 cm⁻¹.
[Adjacent review](larger_RHF16_native_review.json) ·
[All-same spacer review](all_same_pair24_native_review.json) ·
[Alternating spacer review](alternating_pair24_native_review.json).
The fourth double-graft case retains its full Hessian and 19 completed gradient
points after an HDF5 checkpoint-file I/O failure. Its G is blank until full
validation and independent review succeed. The [source-bound saved-matrix recovery](saved_matrix_recovery_01/README.md) retains the full response and requires fresh physical validation. Raw followups are preserved at
`/nfs/lambda_stor_01/homes/rzhu/redox_candidate_continuations/2026-11-02_viologen_larger_RHF16_followups/`.

The article/SI compile at 20/33 pages from the immutable
`2026-11-02_three_ionic_cycles` epoch. The deck has 16 main and 6 supporting
slides: slide 13 compares electronic interaction definitions and slide 20 reports
the three G cycles and their limits. Earlier snapshots and failures below remain
historical evidence. Publication readiness is incomplete.

## Full primary molecular control

[Cycle reconstruction](primary_cohort_cycles.json) · [Scalar/context inputs](primary_cohort_scalar_inputs.json)
· [Source summary](primary_cohort_summary.json) · [Component QC](primary_native_component_reviews.json)

Nine molecular Hessians were freshly computed under the corrected screening
implementation at the exact source geometries. Each own projected nuclear and
ground-spin partition is combined with independently checked corrected
ωB97M-V/def2-TZVPD electronic energy in SMD MeCN. Chloride has an independently
reconstructed monatomic ideal-translation term and no internal modes/Hessian.
The same original r2SCAN-D4/SVP geometry/Hessian protocol and 25/50/100-cm⁻¹
rotor definitions are retained. No geometry repair or ensemble inference is used.

All 81 native SCF checkpoints pass actual molecule, charge/spin, basis and
coordinate-frame checks. Integer electron populations, occupied-orbital
orthonormality and global S² are reconstructed from saved orbitals and AO overlap
matrices. The largest spin contamination is 0.03190, below the original 0.1 gate.
The producer rounds diagnostics to four decimals; comparison respects that
precision, while the independent gate uses full precision. A separate
[70-wavefunction audit](electronic_wavefunction_spins.json) passes for all 29
electronic-control execution namespaces, including their selected parents.

The six complete cycles cover two Fc-referenced reductions, three fixed
precursor-specific quaternizations and disproportionation. Maximum changes from
the retained source across all cutoffs are 2.43963e-8 V and 1.54199e-6 kJ/mol.
The maximum Cartesian-Hessian and nuclear-term changes are 2.22942e-8
Ha/Bohr² and 2.25773e-8 eV. These are measured numerical repeat differences,
not experimental error estimates. No HCl/base term is added to quaternization.
No speciation, kinetic/yield, exact-graft experimental or material accuracy is
inferred, and other benchmark/base-reaction nuclear coverage remains separate.

Raw namespaces and independent reviews remain at:

- `/nfs/lambda_stor_01/homes/rzhu/redox_kernel_controls/2026-11-02_exact_viologen_primary_native_panel_v2/`
- `results/validation/primary_native_spin_checked_2026-11-02/`
- `/nfs/lambda_stor_01/homes/rzhu/redox_kernel_controls/2026-11-02_exact_viologen_primary_native_cohort_checked/`

The first preparation's missing archive hashing import was caught before
scientific execution. The second freezes the shared hashing helper explicitly.
Largest closed/open-shell pilots passed full physical/resource review before
release of the other seven states. The native reviewer also caught and corrected
an import-order error and a comparison against rounded spin metadata; neither
changed a quantum calculation or physical threshold.

## Ionic polymer components and pending cycles

[Current eight-component collection](partial_cycles_08.json) · [Summary](native_collection_summary.json)

Six charge-2 singlet single grafts and two matched neutral empty references pass
independent native review. Lowest internal frequencies range from 7.90 to
12.74 cm⁻¹. The common read-only collector preserves exact graphs, atom/charge
balance, method, solvent/model, thermal definition, ion convention and source/QC
information at every arithmetic boundary. Four whole-chain occupancy cycles
remain blank because complete matched native partners are still unavailable.
These alternate smooth-surface MeCN r2SCAN/SVP values cannot be mixed into the
primary molecular or older primary polymer cycles. No isolated motif binding,
electronic coupling, ensemble or experimental spacing distribution is inferred.

## Larger resource recovery

The original 137-atom charge-4 pilot stopped during RKS H1 construction at an
internal 6-GiB DF allocation limit, with an observed GPU peak of 25,178 MiB.
Its converged SCF, actual grids/orbitals and completed solute partial tensor are
retained. This is an allocation failure, not a failed minimum/derivative gate;
it supplies no native G.

The [verified recovery start](RHF16_recovery_start_verification.json) uses a fresh
namespace, the exact saved frame and partial tensor, a prospectively declared
16-GiB RKS workspace and 40-GiB compact execution budget. The H100-class minimum
and 88,000-MiB observed-peak ceiling remain explicit. The restored-density energy
reproduces exactly (0.0-Ha difference). Only the solute partial tensor is reused:
full solution H1 including PCM, original-equation response, complete Hessian,
all sampled gradients and minimum/nuclear checks remain mandatory. No physical
gate is relaxed. The source failure and prelaunch dependency refusal remain intact.

The [three matching followups](RHF16_followup_preparation.json) are prepared under
the new resource class and remain held for its independently checked completion.
The older prepared class is preserved separately. The current recovery is at
`/nfs/lambda_stor_01/homes/rzhu/redox_candidate_continuations/2026-11-02_viologen_pair23_partial_RHF16_shared/`.

## Workflow, provenance and portability

Shared CPU reviewers now enforce execution closure, exact ownership, finite
resources, actual SCF frames and orbital-spin checks. Their as-used source bytes
are captured and frozen beside receipts. The existing 29-execution/70-frame
panel passes the shared closure/frame review with unchanged electronic values.
A shared bounded dispatcher serves both native workflows with authoritative
reservation discovery, live node preflight, disjoint CPU/thread caps and aggregate
RAM/storage accounting. A fast observation timeout can trigger one longer query
through the same gate; it never restarts a worker or reclaims a live claim.

Quantum preparation binds a declared shared GPU4PySCF source root only when
its grid, integrator and compiled-kernel bytes match the installed environment.
The discovery/config is frozen in the namespace. This resolves node-local prefix
differences without changing the executable, method, libraries or source data.
The current focused suites and compiler/visual/export records are retained;
public data DOI, final reuse terms and author metadata remain outstanding.


## Versioned corrected CPU example

The [new standalone example](../../examples/exact_viologen_corrected_replay/README.md)
contains ten states, nine new Hessians and all nine saved-gradient panels.
The [fresh extraction check](corrected_CPU_validation.json) passes with only
hash-locked NumPy 2.2.6 and six original source/calculation/prepared-package
prefixes blocked. Maximum frequency/nuclear replay errors are 6.49e-11 cm⁻¹
and 3.60e-14 eV. The deterministic 659,642-byte archive has SHA-256
`d4b0f13c9b845c4f42090b9802a94d683be70c111313b2768baa3c8b86c2502e`.
This cache does not recompute SCF, orbital spin or resource qualification;
those have separate checked raw records. The original six-panel example is intact.
A missing compression-report argument was caught after archive creation; the
archive was preserved, its report recovered and every member checked by fresh
extraction. The CLI now rejects missing arguments before creating an archive.
