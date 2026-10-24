# Canonical GPU energy after checked geometry

The [plan](plan.json) combines an independently reviewed V100 resource point
with the exact stationary optimizer output. It computes fresh gas/SMD
ωB97M-V/def2-TZVPD energies on that geometry through the canonical SP writer.
Both source parents and their hashes remain explicit. The resulting record can
be adopted by the native stage without another hybrid calculation.

The first atomic reservation was refused before calculation; its launch and
log remain preserved. The [actual start](startup_verification.json) is verified
on lambda2 GPU 7 with disjoint CPUs 5/19, matching thread caps, a 283 GiB RSS
bound and 32 GiB free-RAM floor. Final energy admission is pending. This stage
performs no optimization, native Hessian or free-energy calculation.

The general CPU/GPU parent join rejects mismatched backend, state, geometry,
protocol and nonfinite or unconverged energy. Tests verify source preservation
and native adoption without recomputation. Raw records and per-SCF recovery
files remain under `calcs/polymer_models/2026-10-24_GPU_electronic/`.
