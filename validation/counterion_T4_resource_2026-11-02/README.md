# Separate T4 fixed-gradient resource control

One exact 42-atom chloride-contact state (+1 singlet) passes an independently
reviewed fixed-geometry primary r2SCAN-D4/def2-SVP/SMD-MeCN SCF/gradient on a T4.
The [review](fixed_gradient_review.json) reconstructs the actual wavefunction,
electron/spin/frame/energy settings and saved gradient, then compares with the
checked V100 reference. Energy and maximum force differences are
1.36e-12 Ha and 1.30e-12 Ha/Bohr. Observed peak GPU use is 9,088 MiB, below the
12-GiB ceiling with at least 25% headroom on the reported 15-GiB device.
The source, environment, runtime FP64 smoke checks, closed ownership and host/GPU
resource chain are retained in the raw namespace and review hashes.

This qualifies one fixed SCF/gradient resource control. It does not qualify
optimization, Hessians, nuclear terms, native G, general T4 molecular work or
material performance. Full molecular/Hessian execution still uses the original
validated 32-GB-class default. A separate optimization pilot is required before
releasing the prepared counterion geometry repairs to smaller cards.

The maintained entry is `python -m redox.workflows.counterion_gradient_pilot
--help`. Device inspection is separate from workload admission: inspecting a
T4 cannot bypass the full-method memory requirement. The [17 scope checks](scope_checks.xml)
reject other molecules/states, optimization/Hessian work and changes to the
strict prospective memory profile. Two allocated CPUs and six native caps=2
remain mandatory. The model, basis, SCF tolerances, grid response and run-local
screening correction match the source reference.

The first adapter attempt failed before GPU/guard/scientific startup because
of a missing plan read in the launcher. Its log and placement remain preserved;
the exact old controller was confirmed closed and science outputs absent before
retry. The second placement completed normally under the standard system gate.
No existing calculation or source result was overwritten. The current header is
[current.json](current.json); bulk and lead Hessian calculations continue with
their established resource classes.
