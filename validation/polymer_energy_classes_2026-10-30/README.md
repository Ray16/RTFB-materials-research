# Remaining occupied-polymer energy workflow

[Audited roster](summary.json) · [Resumable queue inputs](queue_inputs.json) ·
[Checks](tests.xml) · [Preparation manifests](inputs/)

All 160 declared inputs are accounted for: twenty inputs for each combination of
two specified backbone configurations and four phthalimide occupancy graphs
(site 3, site 4, sites 2/3 and sites 2/4). There are 156 prepared electronic
calculations and four independently checked resource points reused only for
their identical scored geometry. The roster checks one chemical graph and one
backbone stereo fingerprint within every group. **Preparation is complete;
quantum evaluation and ensemble thermochemistry are not.**

Six failed geometry attempts remain preserved. Their independently checked
recoveries retain exactly the original declared task and input identity; five
require new energies and one already has its exact checked energy. A previous
ten-state preparation is also reused only after matching its task, scored
coordinates, solvent method, numerical settings, runtime and resource class.
No geometry or numerical result is transferred to another conformer.

The energy jobs retain the archived scientific worker in a read-only source
view. Measured same-constitution resource evidence sets the 160 GiB single-site
or 240 GiB paired-site host-memory guard, with at least 25% headroom over the
measured peak and a 32 GiB free-memory floor. GPU execution requires the proven
32 GB class and the system reservation gate. The per-state launcher allocates
disjoint CPU pairs and all six native thread limits. The dispatcher additionally
accounts for unused RAM promises and unresolved launch handoffs; currently small
RSS does not make that promised memory available for another worker.

The resumable dispatcher uses atomic per-case claims, preserves uncertain RPCs
and partial scientific outputs, and stops a host after three consecutive actual
scientific failures. It does not automatically reopen a failed scientific state.
A separate collector applies the unchanged independent electronic-energy audit
on the source host after the worker and its residency guard close. Completion of
an electronic state never supplies a missing DFT Hessian or nuclear partition.

[Availability discovery and gate-contention recovery](discovery_recovery/README.md)
preserve verified peer nodes when another observation times out. A gate refusal
can be retried only after its exact producers are stopped and the absence of
any scientific stage, guard or raw output is established. Existing scientific
calculations are not reopened by a discovery timeout.

Original geometry ownership is checked by exact host/PID/process-start identity.
The remote query uses system Python and disables SSH connection multiplexing,
so it requires neither a moved interpreter nor a writable control socket inside
the read-only source view. Unknown ownership remains pending. Dependency refresh
retains prior manifests and existing plans and never rewrites a running or
completed namespace. Forty-seven focused checks cover resource headroom, scope,
ownership failures, recovery identity, source binding and queue admission.

The energy-to-native handoff also preserves original plans and independent
reviews. When a new energy plan contains paths absent from the archived source
view, a separate derivative relocates only its source-map path keys. Original
bytes, exact coordinates, chemical fields, electronic values and numerical
settings remain bound and unchanged. The normal native preparer then verifies
the complete relocated parent before constructing a fresh Hessian calculation.
This is provenance handling, not another electronic calculation or a relaxation
of physical acceptance criteria.

Live status is at
`results/validation/energy_classes_2026-10-30/queue/status.json` and
`results/validation/energy_classes_2026-10-30/independent_reviews/status.json`.
The initial explicit placements and detached controller receipts are retained in
the same directory. Raw calculations and source snapshots are outside the
checkout at `/nfs/lambda_stor_01/homes/rzhu/redox_polymer_continuations/`.
The namespace paths and hashes are recorded per state in the roster.

To run a new dispatcher window after verifying that its preceding owner stopped:

```bash
PYTHONPATH=src OMP_NUM_THREADS=1 OMP_THREAD_LIMIT=1 OPENBLAS_NUM_THREADS=1 \
MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 XTB_THREADS=1 \
python scripts/pipeline/run_namespace_energy_queue.py \
  --manifest validation/polymer_energy_classes_2026-10-30/queue_inputs.json \
  --out results/validation/energy_classes_2026-10-30/queue \
  --python /absolute/path/to/the/qualified/python --hours 24 --resume
```

The workflow is a source-preserving continuation on this cluster, not a claim
of unrestricted deployment portability. The separate cached CPU example is
self-contained. Final native free energies, duplicate analysis, complete
10-to-20 comparisons and balanced occupancy ensemble cycles remain required.
This separate smooth-surface acetonitrile control cannot supply replacement
components for the primary-model PMDI or phthalimide cycles. The manuscript
retains its checked 67-state snapshot until new scientific comparisons qualify.
