# Recorded scientific environments

`quantum-linux-64.explicit.txt` contains exact conda package URLs and MD5 package hashes for the Linux environment used for the checked calculations. `quantum-pip.lock.txt` pins all 236 distributions installed through pip. `quantum_runtime.json` records the interpreter, critical versions and qualified libxtb SHA-256. The separate CREST sampling environment is recorded under `results/validation/conformer_revision_2026-10-22/environment_explicit.txt`.

A prospective restoration starts with `conda create -n NEW_ENV --file environments/quantum-linux-64.explicit.txt`, then installs the pip lock using the CUDA 12.8 PyTorch wheel index for the `+cu128` torch build. The ordinary PyPI index is needed for the other packages. CUDA drivers and a supported device are external requirements. These are an explicit conda artifact lock and a pip version lock; pip wheel artifacts are not yet hash-locked. A fresh installation of the complete GPU environment has not been tested, and the primary deployment still contains cluster paths. Do not describe these files as a portable, qualified GPU release.

The bootstrap `setup_env.sh` now pins xTB 6.7.1 and its exact build, and pins the five critical scientific Python libraries to their observed versions. Its default convenience installation is distinct from the recorded explicit environment. The installed xTB library must still match the primary qualification hash; changing or disabling that check does not establish a new qualified environment. Path configuration and corresponding certificate requalification remain release work.

For an immediately portable example requiring only standard Python, use `examples/cached_thermochemistry/reproduce.py`. It reconstructs a potential and two grafting sums from ten cached state records without importing any of the quantum packages or accessing the original paths.


The [scientific review package](../results/releases/primary_review/README.md)
provides a relocated standard-library replay of 243 thermochemical states and
32 electronic-only records, including primary and auxiliary SI quantities. Its
optional CPU check reconstructs 242 Hessian mode sets and nuclear sums, plus
restart geometry/frequency comparisons. It loads no libxtb and opens no original
calculation paths. Explicit environment locks are included; a fresh GPU
installation and portable primary deployment still need qualification.


GPU4PySCF's import name is `gpu4pyscf`, while the installed distributions are
`gpu4pyscf-cuda12x==1.8.1` and `gpu4pyscf-libxc-cuda12x==0.8.1`. Both are already
in the pip lock. A failed metadata lookup for the import name does not imply
an editable or unrecorded installation; runtime checks inside reservations
still verify the imported implementation.
