# Qualified xTB 6.7.1 runtime (Linux x86-64)

This explicit 18-package conda closure pins the thermodynamic kernel and CLI
used by the checked native thermochemistry: `xtb 6.7.1 gfortran_ha02e54d_5`.
The package URLs and MD5 hashes are in `environment.explicit.txt`; full artifact
SHA-256 values, exact dependencies and qualified library/executable SHA-256
values are in `packages.json`. This is the xTB component environment, not a
complete Python/CUDA or molecular-DFT environment.

```bash
conda create --prefix /absolute/path/xtb-runtime --file environment.explicit.txt
export REDOX_XTB_PREFIX=/absolute/path/xtb-runtime
```

Set the prefix before importing thermochemistry modules. `qrrho`,
`native_thermal` and `projected_thermal` now share this configuration. With no
explicit prefix, discovery checks the Python interpreter and base-interpreter
prefixes; the historical cluster prefix is a compatibility fallback. A fresh
installation therefore does not require that historical path.

For a nonstandard layout, `REDOX_XTB_LIBRARY` and `REDOX_XTB_EXECUTABLE` override
the corresponding files individually. Explicit overrides must be absolute,
nonempty paths. The qualification hash remains mandatory for the native kernel:
path relocation is independent of changing the scientific implementation.
`qualified_xtb_prefix` remains available for an explicit process-local binding
before a thermodynamic worker or Python threads have started.

The [fresh-prefix check](../../validation/xtb_runtime_relocation_2026-10-31/README.md)
installs all 18 exact packages in a new temporary directory, checks installed
artifact metadata and both qualified file hashes, and checks dynamic linkage
against that prefix and the system libraries. It runs the actual thermodynamic
fixtures with the original library/executable paths masked. The existing Python
environment is used; complete fresh Python/CUDA deployment remains separate.

To repeat this component check with a conda installation:

```bash
python scripts/release/check_xtb_installation.py \
  --lock environments/qualified_xtb --manager /absolute/path/conda \
  --out /new/durable/check-directory
```

Use `--offline` when the locked packages are cached. The optional `--bwrap`
argument masks the original xTB paths where those original files exist. No GPU
calculation runs. Failed records are preserved; only the newly created temporary
installation is removed. `lock_xtb_runtime.py` can reproduce the lock from an
existing qualified installation without modifying that installation.
