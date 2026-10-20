# Grafting submission and maintained visual results

The finite MeCN campaign contains 83 new component calculations for the 75
balanced local attachment reactions: 72 unique precursors, three shared references
and eight product states. Another 67 product free energies are reused with
source fingerprints. The broad autonomous research goal remains paused; this is
the explicitly requested grafting submission and visualization work.

Ten reserved GPU workers were independently verified: eight H100s on lambda13
and two V100s on lambda11. Verification checks live PID/start identities, gate
ownership, one resident process per assigned GPU, disjoint two-CPU allocations
and native thread caps. No 16-GB card runs this molecular DFT/Hessian workflow.
The [plan and ownership evidence](../../../grafting/campaigns/2026-10-19_mecn_reference/plan.json)
remain distinct from qualified result counts.

Open the [visual report](../../../grafting/index.html). Qualified reactions appear
in a maintained ranking separated by reference-reaction class. Coverage, exact
structures for all 75 candidates and three vector reaction schemes are included.
The finite workers collect completed records and refresh the report after each
item. The [current status](../../../grafting/status.json) owns the changing counts;
the [first-results checkpoint](first_results/current.json) is immutable.

Validation at this milestone:

- Nine targeted software tests pass, including all 75 attachment identities,
  independent atom/charge balances, the oligomer control cycle and reaction-energy signs.
- Eight visual/browser checks pass, covering 75 structure cards, filtering,
  desktop/mobile layout, image paths and vector molecular exports.
- Ten first-checkpoint reaction free energies are independently reconstructed
  from active SP and admitted thermal components, with matching units and signs.
- Five PDF checks pass. The 13 main research slides retain their text; two
  grafting checkpoint frames are added in Backup on pages 15–16. The checked
  deck has 69 pages. The slide checkpoint is dated; the visual report continues
  updating with the queue.

Scientific scope remains explicit: N–H/O–H reactions retain neutral HCl; N
quaternization retains chloride. Base capture is uncomputed. These are local
capped-site reference free energies, without a validated prediction of rates,
yields or actual polymer behavior. Geometry, spin, native Hessian and thermal
failures are preserved; missing quantities remain blank. Atomic chloride uses
the qualified isolated thermodynamic kernel, independently checked against the
monatomic analytic partition, without inventing vibrational modes.

Bulk inputs and outputs remain on NFS under
`calcs/grafting/2026-10-19_mecn_reference/` and the campaign workspace. The frozen
plan binds inputs and numerical producer hashes; per-item completion receipts
bind the actual output records. Git records the code, plans, launch verification,
maintained visuals and this checked initial result milestone.
