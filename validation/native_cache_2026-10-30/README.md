# Portable numerical inputs for the 75-state polymer evidence snapshot

The [isolated archive check](isolation.json) reproduces all 75 states in the
current paper's expanded phthalimide-control snapshot. It creates a fresh
CPython 3.11 environment, installs only hash-locked NumPy 2.2.6, extracts the
archive, and blocks Python access to the current checkout, archived project,
original calculations and prepared cache. No project package, ASE, quantum
program or xTB library is installed in that environment.

The [archive](</nfs/lambda_stor_01/homes/rzhu/redox_release_candidates/2026-10-30_native_inputs75_checked.tar.zst>)
is 63,350,888 bytes and contains 2,256 files. Its SHA-256 is
`19315d6a19d048fe249174ba31549fed72d99ccf141de26a9401c505c8d254e7`.
The [package manifest](package_manifest.json) records content identities and
per-state provenance; [compression.json](compression.json) binds the archive.
Bulk data stay on shared storage. Public archival deposition and final metadata
remain pending.

The cache supplies geometries, Cartesian Hessians, gradients, sampled
displaced-gradient SCF records, electronic parents and three-cutoff thermal
partitions. Its [replay](replay.json) reproduces the saved derivative checks,
rigid-space projections, internal modes and nuclear/free-energy sums. Maximum
differences are 6.97e-11 cm^-1 for internal frequencies and 1.50e-13 eV for nuclear
terms. These are reconstruction errors, not estimates of chemical accuracy.

The unchanged original reviews retain hashes of historical executable binaries
and response scratch tensors omitted from this numerical cache. Complete
historical source closure and live GPU execution are not revalidated. The cache
does not repeat SCF, supply the missing occupied conformers, prove global
sampling convergence, or qualify a full fresh GPU deployment.

The 75 records comprise forty empty states, 18/16 site-2 occupied states across
the two backbone configurations and one site-4 resource point. They match the
existing immutable scientific snapshot; this release introduces no new energies
or candidate-admission decision. Newer completed records remain separate until
their physical and sampling analyses are integrated.

Run the extracted archive with:

```bash
python3.11 -m venv ../native-replay-env
../native-replay-env/bin/python -m pip install --require-hashes -r requirements.lock.txt
../native-replay-env/bin/python -B -I replay.py replay --package . --out replay.json
```

The environment lock targets CPython 3.11 on Linux x86_64. For a dependency-free
potential/grafting arithmetic example, use the separate
[ten-state example](../../examples/cached_thermochemistry/README.md).

Maintainer preparation, packaging and isolated validation are implemented in
[native_cache.py](../../scripts/release/native_cache.py) and
[validate_native_cache_isolation.py](../../scripts/release/validate_native_cache_isolation.py).
The saved [environment log](environment.txt) includes installation and dependency
checks. A small [real-state fixture](../../tests/data/native_cache_one_state.tar.gz)
supports changed-data and rehashed-but-inconsistent-input tests; its provenance
is recorded in [fixture.json](fixture.json). The focused cache and paper tests
pass, as recorded in [tests.xml](tests.xml) and
[integration_tests.xml](integration_tests.xml). The manuscript checker compares
the portable results with the same frozen native states and rejects unsupported
coverage or release claims.
