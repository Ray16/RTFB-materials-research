# Exact-PMDI neighbor-cycle CPU replay

`reproduce.py` is packaged with fourteen original cached states by
`scripts/release/prepare_neighbor_addendum.py`. The prepared archive's location,
checksums and relocated replay are maintained in
[the release index](../../results/releases/neighbor_cycle/README.md).

After extracting the archive, run its copy of `reproduce.py` from any directory:

```
python3 -I -S /path/to/package/reproduce.py --data /path/to/package
```

Only Python's standard library is required. Original absolute source paths are
descriptive; no original project file is opened. The calculation reconstructs
one molecular potential, fourteen balanced reaction comparisons, and the finite
PMDI conditional neighbor cycle at three rotor cutoffs. Raw Hessians and mode
records are included, but this arithmetic example does not reintegrate nuclear
terms, run quantum chemistry or establish conformer convergence, isolated motif
binding, electronic coupling or synthesis eligibility.
