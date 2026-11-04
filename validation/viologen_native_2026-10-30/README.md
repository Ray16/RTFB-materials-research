# Native continuation of an exact enumerated viologen

The baseline enumerated design `search_viologen_b7a4dca7a5e1` is the same exact
capped graft as the registered `viologen` model: InChIKey
`PBNXFHGCSKYNHC-UHFFFAOYSA-N`. Its independently reviewed Merrifield attachment
route passed before this campaign. The campaign computes the missing native
thermochemistry for the declared (+2, singlet), (+1, doublet) and (0, singlet)
states. It does not infer an experimentally feasible material.

All three states have geometry-matched, independently revalidated
wb97m-v/def2-TZVPD/SMD electronic parents. The older optimizer JSON can still
contain legacy thermal or basis labels; the addressed energy record and selected
native thermal policy govern analysis. No legacy GFN2 correction or xTB Hessian
is substituted for a native DFT term.

At this checkpoint, the dication has independently checked native free energies
at 25, 50 and 100 cm−1 rotor cutoffs. The radical-cation and neutral source
geometries fail the native minimum gate, with initial internal frequencies
−51.79 and −12.04 cm−1. Their bounded signed-mode repairs are running. The
original geometries, failed Hessians, subsequent SCF/gradient evaluations and
optimizer-returned coordinates are retained. Complete reduction potentials are
published only after the requested state set and matched reference pass.

The normal scientific implementation remains source-frozen: r2SCAN-D4/SVP/SMD
native derivatives and projected qRRHO with the original convergence, gradient,
Hessian, sampled-column, frequency and spin gates. The repaired geometry's
active single point must match its actual coordinate hash. The independently
reconstructed Fc reference uses the same active/native protocols. Its raw
partitions reproduce all three stored reference free energies exactly, with
78 source bindings checked. The spin-only Fc convention and benchmark's
reference/charge-state limitations remain explicit.

The continuation mounts the original archive read-only. Only private molecular
DFT, lock and repair trees are writable at the original recorded paths. Each
state uses the shared GPU gate, a suitable 32-GB-class V100, two disjoint CPUs,
six thread caps, an explicit contraction runtime and durable NFS outputs. A
prospective 16-GiB host guard retains a 32-GiB free-memory floor; this is a
per-state resource pilot, not blanket qualification of that bound for other
molecules or of 16-GB GPUs.

The checkpoint adapter preserves normal method binding when geometry scanners
copy mean-field objects. It saves SCF settings/densities, returned gradients and
Hessians and optimizer-returned coordinates. Actual restricted and unrestricted
CPU calculations reproduce energies, gradients and Hessians with the adapter;
both geometry scanners evaluate the moved geometry correctly. Forty-one focused
checks also reject changed graph/state/basis, missing/nonfinite data, incomplete
sampled columns, stale namespace caches and invalid redox/reference arithmetic.

The independent collector waits for exact scientific ownership and all execution
records to close. It reconstructs displaced-gradient Hessian checks, projected
nuclear terms, state sums, density spin and resource evidence. Missing or
unverified quantities remain pending. A passed process exit alone supplies no G.

| Artifact | Location |
| --- | --- |
| Frozen input/state snapshot | [inputs](inputs) |
| First checked state | [ox2_independent_review.json](ox2_independent_review.json) |
| Initial minimum failures | [initial_geometry_diagnostics.json](initial_geometry_diagnostics.json) |
| Reference reconstruction | [Fc_reconstruction.json](Fc_reconstruction.json) |
| Live independent outcomes | `results/validation/viologen_native_2026-10-30/independent_reviews/status.json` |
| Raw states/checkpoints/repair history | `/nfs/lambda_stor_01/homes/rzhu/redox_candidate_continuations/2026-10-30_viologen_native/` |

Preparation, preflight and guarded launch are implemented in
`src/redox/workflows/molecular_native_namespace.py`; the numerical worker and
collector are `scripts/pipeline/run_molecular_native.py` and
`scripts/pipeline/observe_molecular_native.py`. The publication utility
`scripts/pipeline/publish_molecular_native.py` requires a complete checked state
set, reconstructs the matched Fc reference, and leaves grafting, occupancy,
coupling and incomplete material criteria blank. Its plotting script is kept
separately in `scripts/plotting/candidates/plot_molecular_redox.py`.

This advances the exact-candidate requirement under M1. It does not complete
actual-condition grafting thermodynamics, polymer occupancy, coupling, global
conformer/spin ordering or the all-axis material decision. The manuscript's
next scientific revision will incorporate complete checked potentials and the
geometry/QC consequences; running jobs are not inserted as paper results.
