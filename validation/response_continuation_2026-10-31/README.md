# Checked response-checkpoint continuation

The full-Hessian workflow now supports raw response-block continuation in a new
namespace. This is an implementation and saved-source audit milestone. No GPU
continuation has been launched from these still-live sources, and molecular
continuation equivalence, full Hessian qualification and native G remain open.

The [read-only source check](source_preflight.json) binds six complete blocks
(288 right-hand-side rows) in each of the V100/H100 controls to their original
orbital state, exact grids, retained geometry and full solvent H1 array. These
are two prefixes of the same 119-atom model, not distinct molecular endpoints.
Both original workers were verified alive at capture; both are ineligible for
continuation. Their running calculations remain untouched. Completed blocks
remain unqualified until the full ordinary physical checks pass.

The workflow requires a proven closed owner, matching PID/start identity and a
closed guard. Only a normally returned full control or an explicitly recorded
wall-budget stop can supply this continuation. An observation timeout, unknown
owner, numerical failure or unresolved shared-storage/output limit cannot reopen
an attempt. Original records are retained in their original locations.

Source binding checks the complete SMD H1 including its PCM contribution,
coordinates, original SCF/orbital frame, grid hashes, contiguous block prefix,
FP64 payloads and solver identity. A restored center resolves back to the
original SCF frame rather than borrowing a new density or different grid.
The new H1 capture also witnesses the actual equilibrium-solvent context.
The previously rejected solute-only/frozen-solvent probe cannot supply blocks.

At runtime, a saved candidate is considered only for matching RHS bytes,
shape, original solver and complete H1/grid frame. The current response operator
is evaluated afresh. Its unchanged absolute residual target is the supplied
solver tolerance divided by 100. A candidate failing that test is discarded
and the original equation is solved. Corrupt backing arrays fail explicitly.
The operator receives a copy so work-space mutation cannot corrupt the candidate
returned to the caller. Trace records preserve the original solver identity
through repeated continuations, and label the checkpoint wrapper separately.

The continued calculation still calls the solvent-wrapped Hessian kernel,
recomputes the PCM H1 term, saves a fresh full H1 and raw response blocks, and
runs the complete matrix, sampled finite-difference, stationarity, projection
and minimum checks. No electronic parent or chemical free energy is assembled
by this control. The independent reviewer cross-checks continuation decisions
against the actual solver records and retains the source arrays/hashes.

A separate prospective eight-hour continuation profile retains the same
80-GiB host, 31,744-MiB observed GPU ceiling, four-GiB durable output budget and
36/32-GiB initial/live shared-storage margins. Original running profiles and
physical thresholds are unchanged. Admission remains through the shared GPU
gate, node preflight and disjoint CPU/RAM allocation.

The [98 focused tests](tests.xml) include changed equations, full-H1/coordinate/grid
mismatches, bad metadata, corrupted payloads, missing prefixes, live/unreachable
owners, reused PIDs, lower precision and in-place operators. A two-generation
CPU continuation reproduces the original solution and original solver identity.
The real array-reader checks pass both saved 119-atom prefixes. These results do
not substitute for a completed molecular continuation/fresh-run comparison.

## Use after the original owner is proven closed

Inspect a source without starting a GPU calculation:

```bash
python scripts/validation/check_response_continuation_source.py \
  --source SOURCE_NAMESPACE --out NEW_AUDIT.json
```

Prepare a separate full-Hessian continuation:

```bash
python -m redox.workflows.screening_control prepare-Hessian \
  --spec ORIGINAL_NATIVE_NAMESPACE --reference CHECKED_NATIVE_REFERENCE \
  --resume-response CLOSED_FULL_HESSIAN_NAMESPACE --out NEW_CONTROL_DIRECTORY
```

Then use the existing node probe and gated launcher under the fleet workflow.
Review normal completion with `review_screening_polymer_hessian.py`. CPU tests,
source-prefix binding and successful submission alone never qualify the final
Hessian. The [manifest](manifest.json) freezes producer snapshots and captured
records; raw molecular calculations remain at the source paths in the audit.
