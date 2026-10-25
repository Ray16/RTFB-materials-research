# Cached Fc and solvent-surface control replay

This small CPU example contains original bytes for six Fc electronic records
and eight small quantum control states (four molecules, two surface models).
Each quantum state includes geometry, active gas/SMD energies, a Cartesian DFT
Hessian/gradient, nuclear data and all 49 baseline/displaced SCF records.

Standard-library electronic/component replay, from any directory:

`python3 -I -S /path/to/example/reproduce.py --electronic-only`

For full cached derivative/mode/nuclear replay, create a Python 3.11 Linux x86_64
venv, install the hash-locked NumPy wheel, then run:

`python -m pip install --require-hashes -r /path/to/example/requirements-cpu.lock.txt`

`python -I /path/to/example/reproduce.py`

The replay checks checksums, contained paths, exact state/geometry/model joins,
Fc electronic gaps/counterfactuals, component curvatures, twelve force components
and two Hessian columns at two displacements per control, and independently
reconstructed modes and nuclear sums at three rotor cutoffs. The comparison
controls impose rotational symmetry number one; this is not the primary-state
symmetry procedure. No original project/calculation path is opened.

This replays stored numerical evidence. It does not run SCF/DFT, load CUDA or
xTB, repeat the physical calculations, validate predictive solvent accuracy,
admit a candidate, establish polymer G or demonstrate ensemble convergence.
Code is MIT-licensed. Final public data-deposition licensing/DOI is pending.
