# Final-CREST input storage inventory

The alternating empty-backbone search produced 4,329 topology/stereo-checked
GFN2/ALPB(MeCN) conformers. Its original queue has twenty declared DFT starts.
These low-level minima are not a DFT ensemble or evidence of ensemble convergence.

The recovered [inventory launch](recovery_01/launch.json) measures the exact
screened FP64 full/range-separated DF arrays for the first start:
55,248,211,656 and 24,059,043,712 bytes. The conservative forecast adds 20% and
24 GiB, totaling 120,938,510,217 bytes. Its observed GPU peak is 12,538 MiB;
no SCF energy, optimized minimum, Hessian or G was computed by the inventory.
Raw inventory records are in
`calcs/validation/conformer_inventory_2026-10-24/empty_alternating_recovery_01/`.

The first probe stopped at its GPU telemetry gate; its original source and
outputs are preserved in the preceding raw directory. The corrected probe
saves shapes and telemetry before testing that gate.

A separate bounded energy/gradient experiment on the same exact start completed
within an 80 GiB RSS bound and 32 GiB free-RAM floor. Its
[independent review](../../reviews/2026-10-24/crest_empty_V100_point/summary.json)
reports a 20,390 MiB GPU peak. This measured resource profile is distinct from
the conservative forecast and does not qualify an optimized conformer or G.
