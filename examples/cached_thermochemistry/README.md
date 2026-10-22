# Portable cached thermochemistry example

Run `python3 -I -S examples/cached_thermochemistry/reproduce.py` from any directory (using the script's absolute path when outside the checkout). Python 3.10 or later and its standard library are sufficient. There are no NumPy, RDKit, PySCF, xTB, CUDA, NFS or credential dependencies. The entire example can be copied to a separate directory or computer.

Ten exact cached states include geometries, active electronic-energy records, DFT-Hessian records, projected frequencies/thermal partitions and a monatomic chloride partition. The manifest binds every file by SHA-256. Original absolute paths are retained only as descriptive provenance and are never opened by the reproducer. The example checks geometry/charge/spin identities, mode counts, cached partition arithmetic and atom/charge balance before evaluating G sums.

It reconstructs the first reduction of the capped PMDI-ethyl model (−1.075790864 V versus the computed Fc couple), its HCl-reference grafting G (−0.291571 kJ/mol), and its DBU separated-ion route G (−103.068826 kJ/mol). New conversions use 96.48533212331002 kJ/mol/eV; the documented historical 96.485 factor is also reported and checked against its source. Both displayed reactions have Δn=0, so their 1-atm-to-1-M correction cancels.

This performs arithmetic from saved, checked energy and thermal components. It does not rerun quantum calculations, reintegrate qRRHO, reproduce GPU execution, independently establish the geometry/Hessian QC or validate chemistry experimentally. Full computational reproduction and a public DOI archive remain separate tasks.

The source export is reproducible with `scripts/release/prepare_cached_thermochemistry.py --out NEW_DIRECTORY` in the original scientific environment; it refuses incomplete backing states and preserves all original bytes. That exporter is not needed to run this example. [Relocation check](../../results/validation/cached_thermochemistry_2026-10-22/summary.json) ran a copied example from `/`, under `/usr/bin/python3 -I -S` without the checkout environment. Negative tests reject changed energy bytes, a rehashed wrong geometry identity and an unbalanced reaction.
