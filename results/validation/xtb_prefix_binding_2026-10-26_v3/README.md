# Explicit xTB prefix binding: scoped relocation check

[Independent review](summary.json) · [Focused tests](tests.xml) ·
[Original-prefix probe](original.json) · [Relocated-prefix probe](fresh.json) ·
[Relocated-interpreter probe](fresh_python.json)

The new `scripts/pipeline/run_portable_thermal.py` entry point accepts an explicit
`--xtb-prefix` or `REDOX_XTB_PREFIX`. It binds the xTB library and executable
process-locally for the existing thermochemistry modules, without editing their
source or replacing their scientific certificates. Both files must match the
previously qualified SHA-256 identities. Binding is refused after a qRRHO helper
has started or when another Python thread is active; the original module paths
and loader environment are restored on normal exit or failure.

Three independent probe processes used the original prefix, the restored prefix
with the original interpreter and the restored prefix with its own Python
interpreter. The kernel's atomic masses and five thermal quantities matched
exactly in this control. Each helper confirmed the actual loaded library path.
A separate reviewer ran both xTB command-line executables on the same previously
qualified formaldehyde geometry and native Hessian. The parsed RRHO contributions
matched exactly. Eleven focused tests cover path escapes, changed bytes, thread
caps, cleanup, worker-state reuse and altered or overpromoted results.

For a CPU-only check, with `PYTHONPATH=src` and a qualified prefix:

```bash
python scripts/pipeline/run_portable_thermal.py --mode probe --xtb-prefix /path/to/prefix
```

Projected mode delegates to the existing projected-thermal CLI; native mode
requires a real `gpu_reserve` owner and an explicitly pinned CPU allocation.
The entry point sets native thread caps before numerical imports. These two
production modes have **not yet been requalified through this new wrapper**.
The test is a Linux path/kernel/CLI control, not a fresh GPU calculation, a new
native state or proof that all cluster paths and certificate dependencies are
portable. Raw probes are retained under
`calcs/validation/2026-10-26_xtb_path_binding_v3/`; the first probe source and
results are preserved separately.
