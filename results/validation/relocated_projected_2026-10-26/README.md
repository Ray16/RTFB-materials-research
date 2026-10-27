# Relocated native-Hessian partition recomputation

[Independent review](summary_verified.json) · [Focused tests](tests.xml)

Two saved, physically checked spin controls (formaldehyde and methyl radical)
were copied with their native Hessians and certificate backing files into a new
checkout outside the original project. The standard native reader refused both
copies because their stored certificate paths named the previous checkout.
The explicit relocation layer permits only that declared source-prefix change;
all certificate hashes and other metadata must still agree.

The new `scripts/pipeline/run_relocated_projected.py` entry point recomputed both
projected partitions using the separately staged, hash-qualified xTB prefix.
All internal frequencies and nuclear free energies at 25/50/100 cm^-1 matched the
original partitions exactly. Source native records and Hessians retained their
original bytes. The original native-reader function is restored when the binding
ends or raises an exception.

Python opens of the original checkout were blocked during recomputation. A
47,499-line file-system trace across the complete process tree records six xTB
CLI executions and use of the explicit qualified library, with zero accesses to
files in the original checkout. Two metadata calls inspected only the original
root directory. Ten adverse-input checks reject altered gradient, stationarity,
spin, derivative, identity, certificate, frequency, software, asymmetry and
undeclared-prefix evidence. Twenty-two focused tests pass.

This checks CPU postprocessing from saved native Hessians. It does not rerun
SCF, geometries or Hessians, requalify the native GPU entry point, demonstrate
experimental accuracy or establish general platform portability.

The source and raw records are retained at:

- `/nfs/lambda_stor_01/homes/rzhu/0_redox_projected_relocation_2026-10-26/`
- `/nfs/lambda_stor_01/homes/rzhu/0_redox_projected_relocation_2026-10-26_strace/`

`staging_manifest.json` / `input_manifest.json` bind every staged source and input
by SHA-256. The checked trace, newly computed partitions, runtime binding and
negative controls remain in the second directory. The independent review records
their exact paths and hashes. The first review is preserved as `summary.json`;
`summary_verified.json` uses the final reviewer with finite-mode validation.

For new isolated checks, use `scripts/release/prepare_projected_relocation.py`
with an explicit source checkout, the two-control raw root and a new destination.
Run the staged `scripts/release/run_relocated_projected_isolated.py` with Python
`-S`, explicit scientific site-packages, the qualified xTB prefix, the declared
old source root and selectors `formaldehyde,methyl_radical`. The command allocates
one CPU and writes to a new output directory. The check is opt-in; ordinary
native readers retain their existing path-sensitive behavior.
