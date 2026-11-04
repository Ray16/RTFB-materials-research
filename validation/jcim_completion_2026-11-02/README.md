# JCIM completion evidence, 2 November 2026

This milestone strengthens the current JCIM working draft and the exact
viologen lead. Submission remains incomplete. The original paper and prior
frozen scientific evidence are preserved in Git and at their recorded paths.

## Interpretation tests

[Five source-bound tests](context_consequences.json) quantify consequences of
transferring properties between contexts while holding the saved calculations
fixed. Four of five complete DBU reaction sign screens differ from their HCl
reference screens, with all ten within-base pairwise orders preserved. Both
PMDI sites have a positive polymer penalty across three cutoffs; their small
site preference changes sign. Material-mass capacity bounds, sampled contacts
at both sequence spacings and asymmetric TEMPO relaxation provide separate
examples. These are not a benchmark of a prior algorithm or experimental
candidate-selection accuracy. Main Table 2 reports the tested transfers.

## Electronic implementation controls

[Reconstructed sensitivity](component_implementation_replay.json) and
[source inputs](component_implementation_inputs.json) bind a ten-state exact
viologen panel and a twenty-one-input base-reaction panel. The panels share two
Fc reference executions, leaving 29 distinct accepted namespaces. The largest
corrected-minus-installed electronic contributions are 1.24e-11 V for the two
viologen reductions, 1.35e-9 kJ/mol for its three local grafting sums and
4.18e-9 kJ/mol across five HCl-reference plus twenty base-assisted sums.
All 70 saved SCF checkpoints pass [actual molecular-frame checks](checkpoint_frame_summary.json).
A failed first placement affected by foreign GPU residency remains preserved;
its unclosed state was not used as an accepted calculation.

These controls recompute E at fixed source geometries. They do not requalify
retained nuclear terms or total reaction G, demonstrate intrinsic DFT accuracy,
or certify broad production use of the corrected implementation. The earlier
71-state benchmark and 119-atom full-Hessian controls retain their own scopes.

## Portable reviewer evidence

The [archive check](component_archive_check.json) verifies every byte of the
330,674,831-byte local control archive, SHA-256
`5af8f5912d03a76c47e12d3b4ce0d4e2bb8d6a3e180adfe127ae3f3b3101c917`.
A CPU-only extraction replay reconstructs 70 SCF frames, selected energies and
all declared electronic cycle shifts with original project/data paths blocked.
It uses NumPy, h5py and PySCF; it does not run SCF, Hessians or GPU kernels or
independently reconstruct wavefunction spin diagnostics. Raw files remain at
`/nfs/lambda_stor_01/homes/rzhu/redox_release_candidates/2026-11-02_component_controls/`
and the sibling `.tar.zst`. The archive is a local review draft: final computed-data
reuse terms, authorship metadata and public deposition remain outstanding.

## Lead calculations and portability

The [charged-native collection](../viologen_native_2026-11-02/README.md) contains
the first independently reconstructed 114-atom +2/singlet viologen native
component, twelve checked stationary components and four matched E cycles.
Seven smaller/matched-empty native followups have exact live-owner observations;
the 137-atom double-graft class has a separate H100 resource pilot. The separate
smooth-surface r2SCAN-D4/SVP MeCN model is retained throughout those components.
No complete occupied G cycle, converged ensemble or material admission follows.

A declared site configuration replaces assumptions about node home paths.
Never-started dependency relocation requires byte-identical replacement files
and proven-closed controllers. New native preparations freeze their review
inputs inside the worker-visible source namespace. Larger-native resource
bounds are shared between preparation, launch, worker and review, retaining
all physical thresholds and leaving the original 32-GiB class unchanged.

## Checks and remaining requirements

The aggregate focused suite passes 131 tests; the additional resource-class
suite is included. Actual node preflight/startup is recorded separately from
physical completion. The 18-page article and 30-page SI build from the immutable
current evidence epoch. Final PDF inspection and self-contained upload records
are linked through the manuscript build/visual/export reports.

The next scientific requirements are the remaining matched native components,
complete balanced free-energy cycles and declared occupied sampling comparisons.
Counterion/speciation sensitivity, exact-graft applicability and accessibility
remain limits on synthesis/material claims. Public release and final author
metadata are required before submission. Current pointers are maintained in
`results/README.md` and `results/INDEX.json`.
