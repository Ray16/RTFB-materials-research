# Isolated V100-compatible tensor runtime

[Runtime definition](runtime.json) · [Two-spin full-Hessian comparison](../../results/validation/cutensor22_full_spin_2026-10-25/comparison.json) · [Runtime guard](../../scripts/fleet/tensor_runtime_worker.py)

This task-local, system-site-packages overlay uses CuPy 13.4.1, cuTENSOR 2.2.0 and fastrlock 0.8.3 with the existing PySCF 2.14.0/GPU4PySCF 1.8.1 installation. It does not modify the primary environment or change molecular methods. The runtime definition deliberately preserves its preparation-time label; qualification is a separate immutable result.

The guard verifies package and library hashes, checks the actual contraction engine, and runs FP64 contractions before entering a scientific worker. It executes inside the reservation/ownership wrapper, in the same process. Library search paths are explicit. Two complete small singlet/doublet controls reproduce the original fallback Hessians to a maximum absolute difference of 7.61×10⁻¹² atomic units. Actual V100 startup and contraction tests pass. This does not certify the unfinished large molecular Hessian.

The installed cuTENSOR 2.7 runtime passed T4 controls but rejected the V100 smoke test with CUTENSOR_STATUS_NOT_SUPPORTED; that negative result is retained separately. A 16 GB fixed-orbital derivative probe is only a subroutine qualification. Full polymer jobs still require the independently qualified 32 GB class and host-memory guards.

Installation records, exact wheels and SHA-256 values remain on shared storage under `calcs/runtime_profiles/cutensor22_2026-10-25/` (`installation.json`, `wheelhouse/`, `venv/`). The overlay is an execution repair with cluster-specific loader paths; it does not close general fresh-checkout portability. No weights or bulk binaries are committed.
