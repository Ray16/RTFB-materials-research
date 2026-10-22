# Repeated real descriptor concurrency

The [checked summary](summary.json) records fifteen actual executions of the same
1,593 recipes / 837 exact designs, with eight tasks per chunk. Each condition ran
three times in a seeded randomized order. Every result was compared with the
independent serial descriptor definition; all 12,555 comparisons matched. Actual
CPU affinities, thread caps, ownership and workers doing useful work were checked.

| Workers | Hosts | Median launch-to-observed completion (s) | Workers doing useful work, min–max |
| ---: | ---: | ---: | ---: |
| 1 | 1 | 4.436 | 1–1 |
| 2 | 1 | 3.147 | 2–2 |
| 4 | 1 | 2.680 | 4–4 |
| 8 | 1 | 2.497 | 7–8 |
| 8 | 2 | 2.585 | 8–8 |

Timing includes node preflight, SSH, scorer initialization, work and up to 0.1 s
of observation polling. Task preparation is measured separately. More workers
help this small library, with diminishing returns: the eight-worker median is
about 1.78 times faster than one worker, rather than an eightfold speedup. These
shared-node observations do not establish an optimal worker count for larger
libraries or quantum work. Host capacities and preflight records remain with
each raw launch. The largest per-worker peak is about 180 MiB; sums of individual
peaks are not simultaneous aggregate memory measurements.

Open the [figure](../../figures/validation/cpu_descriptor_concurrency.pdf). Its
ranges show observed minima and maxima, not confidence intervals. The primary
elapsed measurement uses the controller clock. Separate cross-host worker spans
are retained with explicit uncalibrated clock uncertainty and are not the basis
of the primary timing comparison.

Raw calculation/control files remain at
`calcs/scale_validation/2026-10-22_concurrency/` (about 35 MB). Their immutable
plans and output receipts identify every exact input and result. The checked
handler, launcher and experiment driver are retained in [source_snapshot](source_snapshot/).
This validates a descriptor-stage concurrency experiment, not million-molecule
chemical throughput, quantum scaling or experimental accuracy.

Reproduce with all six native thread caps exported as one:

```bash
PYTHONPATH=src python scripts/validation/benchmark_cpu_concurrency.py \
  --config config/cpu_concurrency_benchmark.json \
  --raw calcs/scale_validation/NEW_CONCURRENCY --out results/validation/NEW_CONCURRENCY
```

The config specifies benchmark conditions, not a global CPU or GPU worker cap.
