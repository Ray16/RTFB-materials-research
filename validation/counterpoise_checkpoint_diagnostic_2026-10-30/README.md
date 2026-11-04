# Preserved AB orbital-stationarity diagnostic

[Startup source bindings](startup.json) · [Launch](launch.json) · [Scope](summary.json)

The previous tight AB SCF refinement stopped with a preserved checkpoint after
63 recorded cycles; it never supplied a completed, accepted term. Its original
producer is confirmed gone. A CPU-only diagnostic now reconstructs the energy
and orbital gradient from that saved checkpoint using the unchanged common grid,
full-complex RI auxiliary basis, frozen nuclei and archived evaluator. The
source archive remains read-only. Eight disjoint CPU cores, 32 GiB sampled RSS
and a 32 GiB free-RAM floor bound the diagnostic; live native-worker RAM promises
are accounted for before starting it.

The computation uses no GPU and does not rerun SCF. Its inputs are rehashed after
evaluation. An eventual small orbital gradient would establish only the tested
saved density's stationarity, not the original SCF convergence, a completed
counterpoise cycle, a geometry minimum, MeCN correction or association free
energy. The full gas counterpoise remains unqualified until all five terms and
their normal independent reviews pass. The 1e-7 orbital-gradient acceptance
threshold is unchanged.

The entry point is
`scripts/validation/diagnose_archived_counterpoise_checkpoint.py`. Its four
pre-computation guard tests are included in the 19-check follow-up report under
`validation/native_followup_collection_2026-10-30/tests.xml`. This milestone
records a running diagnostic; inspect the exact owner and live status before
claiming completion or submitting another attempt.
