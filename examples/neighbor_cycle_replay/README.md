# Exact-PMDI neighbor-cycle CPU replay

`reproduce.py` is packaged with fourteen original cached states by
`scripts/release/prepare_neighbor_addendum.py`. The prepared archive's location,
checksums and relocated replay are maintained in
[the release index](../../results/releases/neighbor_cycle/README.md).

After extracting the archive, run its copy of `reproduce.py` from any directory:

```
python3 -I -S /path/to/package/reproduce.py --data /path/to/package
```

The default mode requires only Python's standard library. Original absolute source paths are
descriptive; no original project file is opened. The calculation reconstructs
one molecular potential, fourteen balanced reaction comparisons, and the finite
PMDI conditional neighbor cycle at three rotor cutoffs. Raw Hessians and mode
records are included, but this arithmetic example does not reintegrate nuclear
terms, run quantum chemistry or establish conformer convergence, isolated motif
binding, electronic coupling or synthesis eligibility.

The current archive also supports independent Hessian/mode and nuclear replay.
In a fresh Python 3.11 Linux x86_64 venv, install the hash-locked NumPy wheel from
the packaged `requirements-cpu.lock.txt`, then run:

```
python -m pip install --require-hashes -r /path/to/package/requirements-cpu.lock.txt
python -B -I /path/to/package/reproduce.py --data /path/to/package --full-native
```

This mode reconstructs thirteen molecular Hessian projections and their nuclear
terms at all three cutoffs, plus monatomic chloride translation. It loads no
quantum program or xTB library and does not repeat the original SCF, force, spin
or displaced-gradient checks. Both modes have been tested after extraction with
original repository reads blocked.
