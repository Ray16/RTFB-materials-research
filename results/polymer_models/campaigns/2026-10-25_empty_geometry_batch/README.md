# Broader empty-backbone DFT refinement

[Declared batch](batch.json) · [Latest dispatch snapshot](dispatch.json) · [New-dispatch hold](dispatch_hold.json) · [Energy dependencies](../2026-10-25_empty_electronic_batch/batch.json)

Thirty-eight further CREST-derived starts cover the empty five-repeat model at two declared backbone configurations. The initial conformer of each configuration is handled separately. All jobs use the existing molecular method and independently check final topology/stereochemistry and total-energy forces. Measured resource equivalence transfers placement bounds, never energies, geometries or certificates. Durable raw geometry evaluations remain in `calcs/polymer_models/2026-10-25_empty_geometry_batch/`.

New geometry dispatch is paused after repeated final force/energy inconsistencies. Existing scientific workers continue and passing geometries advance to canonical gas/SMD energy evaluation through a separate dependency-gated batch. Failed conformers remain excluded and their checkpoints preserved. Foreign GPU attachments stop only the owned worker; reviewed continuations occupy new campaign directories. The dispatch JSON freezes when the coordinator stops; individual stage records are authoritative for workers that finish later.

The [CPU decomposition](../../../validation/conformer_force_diagnostics_2026-10-25/conf005_CDS_decomposition.json) isolates a representative mismatch to the legacy SMD CDS term's derivative-scale behavior. A separately declared same-energy geometry-repair pilot uses broader search directions but retains the original independent final force tests and adds all-axis CDS consistency. It must pass before repair fan-out. The first pilot was interrupted by foreign GPU residency; the [reviewed retry](../2026-10-25_CDS_geometry_repair_interruption/batch.json) preserves that attempt. No tolerance was relaxed.

The [first-site local gradient review](../../../validation/conformer_force_diagnostics_2026-10-25/first_site_local_gradient_review.json) finds agreement between legacy and local energy derivatives for the earlier empty, capped and two PMDI single-site states, while retaining larger-step sensitivity. It does not establish global numerical or conformational convergence.

No new native free energies, converged conformer partition or spacing ranking follow from these geometry/energy stages. The twelve inherited phthalimide native states remain a separate finite two-pose baseline.
