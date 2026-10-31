# Polymer recovery and response-execution controls

The [working census](native_census.json) now contains **121/240 independently
reconstructed native input states**. Six additional states verify 14,643 source
bindings in the [new replay](../../results/validation/response_diagnostic_2026-10-31/native_replay.json).
The maximum sampled derivative relative error is 5.506e-6 (0.000551%) and the
lowest internal frequency is +10.23 cm^-1. This includes the successful
continuation of the previously interrupted `all_same__occ_4__conf005` state.
All six use the same explicitly identified smooth-CDS MeCN control model.

| Configuration | Empty | Site 2 | Site 3 | Site 4 | Sites 2,3 | Sites 2,4 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| All same | 20/20 | 19/20 | 14/20 | 13/20 | 7/20 | 10/20 |
| Alternating | 20/20 | 17/20 | 1/20 | 0/20 | 0/20 | 0/20 |

The [coverage figure](../../results/figures/polymer_models/2026-10-31_response_recovery/native_occupancy_coverage.pdf)
shows these twelve groups. Input counts include duplicate minima and do not
establish complete occupied ensembles, balanced ensemble cycles or material
eligibility. The article retains its separate 75-input numerical comparison.
M1, M4, full M7 and final M9 remain open.

## Allocation stall and a distinct diagnostic fault

The four older site-2 jobs remain live. Their accumulated CPU time is dominated
by kernel time despite abundant globally available RAM. An owned read-only
attach/backtrace/detach left the old worker running and untraced. Although its
raw stack lacked symbols through the namespace, mapping its top instruction
against that same live process's executable mappings identifies NumPy's
`DOUBLE_add_AVX2`. See [old_worker_PC_mapping.json](old_worker_PC_mapping.json).
This locates one observed active instruction; it does not establish every
worker's hot function or a unique causal mechanism by itself.

The reusable [allocation comparison](../../scripts/validation/compare_numpy_allocation.py)
runs fresh one-CPU processes off/on/off against identical stored data. With
`NUMPY_MADVISE_HUGEPAGE=0`, the orbital read and finite-value scan complete in
0.209/0.209 s with the exact source orbital checksum. With advice enabled, the
20-s observation ends before completion. The
[comparison](../../results/validation/response_diagnostic_2026-10-31/allocation_comparison/comparison.json)
retains the censored sample, settings, stacks and hashes. NumPy documents this
import-time option; the kernel documents synchronous reclaim/compaction for
advised huge-page regions. [NumPy](https://numpy.org/doc/2.0/reference/global_state.html),
[kernel](https://www.kernel.org/doc/html/v6.7/admin-guide/mm/transhuge.html).
These host-specific controls support an allocation-stall explanation without
changing machine-wide policy or quantum-model parameters.

Three new traced V100 response probes reproduced the fixed-density source energy
to stored precision, then terminated with fatal signals. The last probe completed
one response block with a true residual norm of 9.10e-8 against its original
1.60e-7 target; the full response did not finish. Those attempts supply no native
G or execution qualification. Their periodic stack logs end mid-frame, and a
fault-address hypothesis is consistent with Python's traceback routine. The
[assessment](instrumentation_assessment.json) states the missing image-base
evidence explicitly. A relevant [CPython issue](https://github.com/python/cpython/issues/116008)
documents crashes while walking changing thread/frame state without the GIL.
Both short CPU stress controls completed under the debugger, so they do not
independently reproduce the GPU-process crash. No numerical-engine fault is
established by this evidence.

The updated helper replaces the asynchronous C watchdog with a normal Python
thread calling `dump_traceback` while holding the GIL. A regression check rejects
any reintroduction of `dump_traceback_later`. Event timing leaves operator inputs
and results unchanged; checkpointing adds explicit device-to-host transfers after
completed solves. The original numerical solver, model, basis and source orbital
frame remain bound by hashes. The normal native memory-policy context is also
applied; a 3-GiB available-memory ceiling reduces automatic response batching.

## Checkpointed control in progress

[checkpointed_start.json](checkpointed_start.json) verifies the revised control's
actual V100 process, reservation, two-CPU binding, six native thread caps and
import-time NumPy option. Its source namespace is
`/nfs/lambda_stor_01/homes/rzhu/redox_polymer_continuations/2026-10-31_response_checkpointed_v4b/`.
The changing `case/status.json`, `case/trace/latest.json`, response-block records
and `placement_40GiB/completion.json` there determine subsequent progress;
startup verification is not completion.

The explicitly selected long profile has a two-hour bound, 2-GiB durable-output
ceiling, 34-GiB initial free-storage requirement and a stricter prospective
40-GiB host-RSS guard. It reuses the existing large inputs and preserves the
production queue's separate 64-GiB gate. Each solved Krylov block is written
atomically with its RHS/state identity and hashes. These are raw solver solutions
before occupied-block reconstruction. Final `mo1`/`e1` arrays and recreated grids
are saved separately if the full response returns. Checkpoints remain unqualified;
the original grid arrays are absent from the old source, and no native Hessian,
partition or material decision follows from this control.

The [latest storage observation](../../results/validation/response_diagnostic_2026-10-31/native_storage.json)
estimates a full retained footprint of 99.91 GiB for the 119 unfinished native
inputs at observed class maxima. Some partial files are already allocated, so
this is neither an additional-write forecast nor a peak bound. Upstream energy
outputs and other users' writes are excluded. Further shared capacity is still
needed for the full campaign.

## Reproduction and scope

```bash
export NUMPY_MADVISE_HUGEPAGE=0
# Use the redox interpreter and cap native CPU thread pools before CPU tools.
python scripts/validation/compare_numpy_allocation.py --raw /path/to/source/raw \
  --out /path/to/new/allocation-comparison

python -m redox.workflows.response_diagnostic prepare \
  --spec /path/to/source/namespace.json --out /path/to/new/control \
  --checkpoint-response
python -m redox.workflows.response_diagnostic probe \
  --spec /path/to/new/control/namespace.json
python -m redox.workflows.response_diagnostic launch \
  --spec /path/to/new/control/namespace.json --gpu AVAILABLE_GPU \
  --out /path/to/new/control/placement
```

Launches use `gpu_reserve`, fresh ownership records, disjoint CPU allocation and
the exact source/runtime preflight. Prior attempts are preserved. The optional
debugger uses an isolated, unprivileged Ubuntu package extraction and a
PID-specific cooperative handshake; shared ptrace policy is unchanged. Its
CPU-only intentional-fault controls test the handshake outside and inside the
read-only source namespace. Debugger binaries are not project dependencies or
redistributed scientific data; [debugger_runtime.json](debugger_runtime.json)
records their versions and hashes. Missing modules were caught in the namespace
CPU control before another GPU launch. The earlier absolute sibling-ptrace
restriction inference was too broad for these owned nested namespaces and is
superseded by the successful attach observations.

Forty-one response/ownership/census/ensemble tests pass, and the new Python files
pass Pyflakes 3.4.0. These checks establish their stated behavior and input guards;
they do not close occupied sampling or qualify all hardware/numerical paths.
