# Corrected screening checkpoint: 2026-10-19

The geometry driver now uses full grid-response gradients and records the
optimizer and SCF convergence results. Native thermal partitions exclude rigid
motions while retaining genuine soft modes, using an isolated installed xTB
kernel. Both primary readers require qualified corrected partitions.

- Full software suite: 391 passing tests.
- Actual CPU geometry-driver smoke: passed; GPU driver qualification is separate.
- Corrected candidate export: 22/25 additional potential pairs; 555 arithmetic
  checks, zero discrepancies, inputs unchanged during collection.
- Corrected global assembly: 28,967 checks, zero discrepancies, unchanged inputs
  during review and a second current-inventory check before publication.
- Matched parent benchmark: seven potentials, four molecules, three studies;
  RMSE 0.1714 V, MAE 0.1471 V. Exact-graft accuracy remains unvalidated.
- Slides: checked 38-page PDF with corrected results and a dedicated new-lead frame.

Remaining work includes three incomplete additional candidates, independent
promotion review for two private NDI minima, solution reorganization and kinetics,
broader experimental anchors, synthesis, solubility/swelling and durability.
Potential-window capacity is bookkeeping, not a validated cycling capacity.

Git contains code, curated inputs, small checked exports and provenance metadata.
Raw Hessians, execution scratch, downloaded datasets and earlier scientific
artifacts remain on shared storage at
`/nfs/lambda_stor_01/homes/rzhu/0_redox/`, addressed by the committed inventories
and certificates. Reproduction of physical admission requires those backing data
and the recorded validated environment. Historical unprojected results remain
preserved and are superseded. Pre-existing staged archive moves are excluded from
this checkpoint and retained in the user's index.
