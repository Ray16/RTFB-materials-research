# Checked polymer continuation states

[Summary](summary.json) · [Physical reconstruction](native_replay.json) ·
[Exact inputs](inputs.json) · [Coverage figure](sampling.pdf) · [Tests](tests.xml)

Seven new neutral phthalimide-control states passed the unchanged independent
review and a second, CPU reconstruction of the native DFT Hessians, sampled
displaced gradients, projected frequencies and three qRRHO partitions. The
reconstruction verifies 10,937 source bindings and retains the original source
hashes through explicit path aliases. No new quantum calculation runs in this
reconstruction. The raw independent reviews, SCF evaluations and matrices remain
on shared storage at the locations bound by each input namespace.

With the sixty retained states, coverage is now 20/20 and 16/20 empty-backbone
states, and 18/20 and 13/20 site-2 states, for the all-same and alternating
backbone configurations respectively. These are states of the separate smooth
acetonitrile model, not replacement components for the primary occupancy cycles.
Counts do not prove distinct minima or ensemble convergence. The missing states,
duplicate review, complete 10-to-20 tests and other occupied graphs remain open.

Five additional electronic-energy states also passed their independent reviews:
two alternating empty backbones, all-same site 3 and sites 2/4 resource points,
and the missing all-same site-2 input 2. Their electronic values, resource peaks
and exact bindings are in the summary. Two previously completed site-4 and
sites-2/3 resource points were separately checked and staged without changing
their bytes. No electronic-only state supplies missing thermochemistry.

The maintained handoff uses a fresh output directory and relocates only paths
in a derived electronic-parent receipt. Original reports, tasks, method settings
and numerical values remain unchanged. Read-only namespaces keep the archived
worker and normal physical gates available after repository refactoring. Every
GPU launch retains the reservation gate, disjoint CPU pairs, six explicit
thread caps, an 80 GiB native host-memory bound and a 32 GiB free-memory floor.
The two larger electronic pilot classes use their declared 160/240 GiB bounds;
the site-2 successor additionally rechecks its measured same-constitution
resource witness and 25% memory headroom. Thirty focused checks cover these
handoffs, staging, source checks and resource-scope restrictions.

To reconstruct these seven finished states without a GPU:

```bash
PYTHONPATH=src OMP_NUM_THREADS=1 OMP_THREAD_LIMIT=1 OPENBLAS_NUM_THREADS=1 \
MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 XTB_THREADS=1 \
python scripts/validation/replay_namespace_native.py \
  --manifest validation/polymer_continuations_checked_2026-10-30/inputs.json \
  --out NEW_REPLAY.json
```

This full physical replay needs the bound raw records and the scientific Python
environment. The separate `examples/cached_thermochemistry` example needs only
the Python standard library and ships its own cached states. Public deposition,
broader polymer sampling and a complete material-selection demonstration remain
publication requirements.
