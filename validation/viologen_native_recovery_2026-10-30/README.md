# Checked redox ladder for the exact enumerated viologen

All three declared molecular states now pass independent native, density/spin,
derivative and resource checks. The exact 4-methylbenzyl capped viologen,
`PBNXFHGCSKYNHC-UHFFFAOYSA-N`, has computed MeCN reduction potentials
**−0.608 and −1.208 V versus Fc/Fc+** at 298.15 K and the primary 50 cm−1 rotor
cutoff. Their spacing is 0.600 V. The corresponding modeled monomer
disproportionation, `2 V+ -> V2+ + V0`, costs +57.91 kJ/mol.

These are unfitted, declared-spin, capped-model results. They do not validate
an exact-graft experimental potential, a globally converged conformer/spin
ensemble, polymer interactions or material feasibility. The all-axis decision
remains pending. Attachment free energy, occupancy and electronic coupling remain
blank in [candidate.json](candidate.json).

[numeric_inputs.json](numeric_inputs.json) contains the five molecular/reference
state sums; [native_QC.json](native_QC.json) records actual native checks.
The radical-cation and neutral source structures had native imaginary modes of
−51.79 and −12.04 cm−1. Their accepted repaired minima have lowest frequencies
+19.94 and +16.02 cm−1. The dication gives +23.73 cm−1. All three use the same
active diffuse-basis protocol and native geometry/Hessian definition as the
matched Fc reference. Nuclear terms come from projected DFT modes, with xTB
used only as the qRRHO integrator.

The original neutral and radical quantum runs completed both signed-mode repair
branches, but canonical promotion encountered `EXDEV`: separate bind mounts
prevented atomic rename despite matching filesystem device IDs. The pipeline
now maps promotion arguments into the single writable case mount. A real CPU
namespace test reproduces the failure and verifies successful atomic promotion
with the mapped paths; source hashes, charge and multiplicity are retained.

Recovery preserves the original failed runs and both verified branches. It
selects the lowest active-energy verified branch after both attempts terminate,
reuses its optimized coordinates and matching single point, and performs fresh
guarded native derivative checks. The neutral fresh run reproduces the saved
verified free energies within 4.2e−9 eV across all three rotor settings
([reproduction](native_reproduction.json)). No expensive geometry optimization
was restarted for this recovery.

The independent density audit also exposed a precision-contract error: the
legacy producer rounds S² to four decimals. The radical's independently
reconstructed 0.7535607357 agrees with the stored 0.7536 at that precision.
The comparison now honors the stored precision while retaining the full-precision
0.1 physical contamination limit. Nonfinite, discrepant and physically excessive
spin values remain rejected. Earlier failed review outcomes remain at their
original paths.

Maintained QC helper sources are copied alongside each review, so they cannot be
mistaken for files inside the older read-only archive. The final collection
explicitly selects the original checked dication and the newly checked repaired
radical/neutral namespaces. Older failures are not relabeled. The guarded
finalizer assembles potentials and numerical inputs only after every selected
state passes; its two-potential [figure](figure.pdf) has been visually inspected.
Bars describe 25–100 cm−1 rotor sensitivity, not statistical or prediction intervals.

The maintained preparation/recovery, numerical worker, collector and publisher
are versioned under `src/redox/workflows/molecular_recovery.py`,
`scripts/pipeline/run_molecular_native.py`,
`scripts/pipeline/collect_molecular_native.py` and
`scripts/pipeline/publish_molecular_native.py`. The finite finalizer is
`scripts/pipeline/finalize_molecular_native.py`; plotting remains separate in
`scripts/plotting/candidates/plot_molecular_redox.py`.

Raw optimized states, SCF checkpoints, Hessians and unsuccessful execution history
remain at `/nfs/lambda_stor_01/homes/rzhu/redox_candidate_continuations/2026-10-30_viologen_native/`
and `2026-10-30_viologen_native_recovery_checked/` under the same parent.
The stable candidate pointer is `results/candidates/verified_cases/viologen/current.json`.
The original [continuation snapshot](../viologen_native_2026-10-30/README.md)
remains a historical checkpoint. M1 still requires attachment/polymer-context
completion; the article and slides await incorporation of this finished ladder.
