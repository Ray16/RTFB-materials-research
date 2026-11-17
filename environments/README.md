# Recorded scientific environments

`quantum-linux-64.explicit.txt` contains exact conda package URLs and MD5 package hashes for the Linux environment used for the checked calculations. `quantum-pip.lock.txt` pins all 236 distributions installed through pip. `quantum_runtime.json` records the interpreter, critical versions and qualified libxtb SHA-256. The separate CREST sampling environment is recorded under `results/validation/conformer_revision_2026-10-22/environment_explicit.txt`.

Use `scripts/release/restore_quantum_environment.py` to create a fresh prefix and an exact archived project checkout. It installs `redox==0.1.0` from that checkout, rather than obtaining an unrelated namesake from PyPI. Torch uses its CUDA-specific wheel index; other third-party distributions use PyPI. CUDA drivers and a supported device remain external requirements.

The exact recorded installation was attempted and preserved. All 235 external distributions installed, but `pip check` found fourteen CUDA-version conflicts with the unchanged Torch wheel. The recorded lock and the running production environment are preserved; the exact replay is not a qualified release. `quantum-pip.cuda128-compatible.lock.txt` is a distinct tested profile, changing only those fourteen CUDA packages to Torch's declared requirements. Its fresh installation passes `pip check`, and both complete spin controls pass the independent quantum comparison; see [the validation guide](../results/validation/gpu_environment_restore_2026-10-23/README.md).

The bootstrap `setup_env.sh` now pins xTB 6.7.1 and its exact build, and pins the five critical scientific Python libraries to their observed versions. Its default convenience installation is distinct from the recorded explicit environment. The installed xTB library must still match the primary qualification hash; changing or disabling that check does not establish a new qualified environment. An [explicit xTB-prefix entry point](../results/validation/xtb_prefix_binding_2026-10-26_v3/README.md) now passes scoped kernel and CLI relocation checks with the qualified binary and library hashes. Requalification of native/projected production entry points and broader cluster-independent path configuration remain release work.

For an immediately portable example requiring only standard Python, use `examples/cached_thermochemistry/reproduce.py`. It reconstructs a potential and two grafting sums from ten cached state records without importing any of the quantum packages or accessing the original paths.


The [scientific review package](../results/releases/primary_review/README.md)
provides a relocated standard-library replay of 243 thermochemical states and
32 electronic-only records, including primary and auxiliary SI quantities. Its
optional CPU check reconstructs 242 Hessian mode sets and nuclear sums, plus
restart geometry/frequency comparisons. It loads no libxtb and opens no original
calculation paths. Explicit environment locks are included; a separate fresh GPU
profile is now checked for two rigid controls; general deployment portability remains open.


GPU4PySCF's import name is `gpu4pyscf`, while the installed distributions are
`gpu4pyscf-cuda12x==1.8.1` and `gpu4pyscf-libxc-cuda12x==0.8.1`. Both are already
in the pip lock. A failed metadata lookup for the import name does not imply
an editable or unrecorded installation; runtime checks inside reservations
still verify the imported implementation.

`quantum-cuda128-artifacts.lock.txt` pins the 235 retrieved external artifacts by direct URL and SHA-256. Install the local project separately from its recorded source checkout. The archived source also needs every hash-required certificate backing file; `scripts/release/certificate_closure.py` validates this closure before calculation. Explicit process-local xTB path bindings retain the qualified hash. The public raw-state release and cluster-independent configuration remain separate from these tests.

The [relocated projected-partition check](../results/validation/relocated_projected_2026-10-26/README.md) now recomputes both saved spin-control partitions through an explicit source-prefix and xTB-prefix entry point. All three rotor-cutoff nuclear terms and internal modes match exactly. A process-tree file trace shows no original-checkout file access, and ten altered physical inputs are rejected. This qualifies the tested CPU postprocessing path; it does not requalify native GPU Hessian execution or every project path.
