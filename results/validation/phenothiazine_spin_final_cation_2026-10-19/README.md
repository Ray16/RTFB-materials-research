# Last experimental-parent cation spin comparison

The plan contains the now-ready DMeOMEEPT radical-cation quartet comparison. Its
previous molecule spin writer finished before preparation. A suitable free H100
NVL on lambda13 passed dependency, host-memory/load and CPU-affinity preflight.
The actual resident worker was verified under the GPU reservation, with two
native threads on disjoint CPUs. GPU memory is not pooled.

The SSH launch observation timed out, but durable submission records and the
subsequent live process/GPU probe confirmed the original worker had started.
No duplicate or replacement calculation was submitted. All outputs and source
hashes remain under this immutable run. `parents/<id>/status.json` records outcome;
startup is separate from physical and benchmark admission.
