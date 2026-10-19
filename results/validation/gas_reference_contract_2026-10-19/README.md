# Strict gas-reference solver contract

Three reference workers now share `redox.qm.gas_reference`: D3TaLES-matched gas,
inner-acene literature, and the generic literature worker. New computations require
finite converged SCF energies, explicit optimizer and final-SCF convergence, and
full quadrature-response gradients. Functionals, bases and fitting settings remain
explicit in each worker. The acene record now reports its actual requested basis.

GPU HOMO arrays are transferred explicitly to host. A missing or positive vertical
anion HOMO withholds the usable electron reorganization value; a diagnostic value
is kept separately with binding QC. An unknown HOMO no longer becomes `False` for
`anion_unbound`.

Cache reuse requires matching molecular/method/reference inputs, producer and
solver hashes, software versions, complete finite outputs, and an `ok` result.
Unqualified, partial or incompatible files are preserved and refused before new
computation. Choose a fresh `--out` (D3TaLES) or `--tag` (literature) for a separately
planned replay. Case writer locks and atomic JSON publication prevent duplicate or
partial output. Failed/incomplete runs cannot emit a completed marker or exit as
successful.

The 514-test suite and actual CPU H2/PBE/STO-3G control pass. The latter
also proves a one-step-stopped optimizer is refused. These are solver controls,
not experimental accuracy, ground-spin, Hessian-minimum or chemical-integrity proof.
An actual CPU water/cation four-point run also passes: the saved Hartree points
reconstruct 0.2900554 eV (reported 0.2901), and a second call reuses the qualified
record without recomputation. See `cpu_worker_smoke/report.json`.
The 482-record read-only inventory found no complete new solver-verification
contract in the saved legacy cohorts. It does not prove their values are wrong and
does not retrospectively certify them. No reference dataset was overwritten or
bulk quantum replay submitted. New records retain the four raw Hartree energies
and optimized geometries so the reported eV difference can be audited separately.

The public `properties/nelsen.py` forwarding change remains pending because that
file is bound by active spin plans. Its only production consumer imports `embed`.
The running spin producers remain byte-identical; the new gas module and edited
reference scripts are outside their bound execution paths.
