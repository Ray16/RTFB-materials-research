# Fresh qualified xTB runtime and path relocation

A fresh temporary conda installation of the exact 18-package xTB runtime passes
real thermodynamic tests with the original xTB library and executable paths
masked. This addresses component installation and automatic path configuration;
it does not close full Python/CUDA deployment portability or public deposition.
No GPU job runs, and the existing Python environment supplies the checker.

The [installation summary](summary.json) verifies the artifact versions/builds
and SHA-256 metadata against the explicit lock, and the installed library/CLI
bytes against their existing qualification SHA-256 values. Dynamic linkage uses
the new prefix and system libraries. The original qualified build files remain
unchanged. The temporary installation is removed after capture.

The upstream thermodynamic fixture reproduces G/H/TS to its 1e-12-Ha limits.
The exact 0.391 cm-1 internal stretch survives geometric rigid-mode projection
and remains in the nuclear partition calculation. The separate thermodynamic
worker uses one-thread native pools, and no CUDA module is imported. These
checks test the actual qualified kernel, not cached arithmetic or a substitute
quantum method.

The maintained qrrho, native and projected entry points now discover xTB from
an explicit `REDOX_XTB_PREFIX` or interpreter-prefix installation. Individual
library/executable overrides support nonstandard layouts; relative or empty
overrides fail. The legacy cluster prefix remains a compatibility fallback.
Numerical constants, model definitions, the qualification hash and physical
admission thresholds are unchanged. All running source-frozen scientific jobs
retain their original code and settings.

The [runtime guide](../../environments/qualified_xtb/README.md),
[explicit lock](environment.explicit.txt), [package inventory](packages.json)
and [installed package record](installed_packages.json) make the component
installation reviewable. The [117 focused tests](tests.xml) cover path binding,
invalid overrides, transitive package closure, malformed/private artifact URLs,
thermodynamic soft modes, unchanged native gates and Hessian comparisons.
[Producer snapshots and file hashes](manifest.json) freeze the tested sources.
The raw check remains at the path recorded by the captured source manifest.

A standalone cached numerical example and its fresh NumPy-only replay remain
separate evidence. This check does not claim a complete fresh quantum stack,
resolve polymer sampling or qualify newly computed reaction free energies.
