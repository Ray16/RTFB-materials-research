# Cata-ADI neutral CDS derivative diagnosis

The first native attempt failed its unchanged Hessian symmetry gate:
maximum asymmetry 5.9612834e-6 Ha/Bohr² exceeds 1e-6. The selected geometry,
active energy and failed Hessian remain preserved; the native correction is
withheld. The [frozen receipt](frozen_receipt.json) verifies the actual worker and
controller are terminal before any derivative replay is prepared.

A host-only evaluation of the installed SMD CDS interface reproduces the full
asymmetry. Subtracting the coarse CDS matrix leaves maximum antisymmetry
7.75e-14 Ha/Bohr². Reducing its finite-difference step from 1e-4 to 1e-6 Bohr gives
CDS asymmetry 7.13e-12. The [stored-matrix review](matrix_review.json) independently
recomputes those quantities and binds all arrays by hash. The
[diagnostic report](cpu_cds_diagnosis/report.json) retains the two worst contact
columns and a derivative-step sweep. No SCF, GPU calculation or canonical
admission is inferred from that host-only diagnosis.

This motivates one [private full finer-step replay](finer_replay/plan.json), using
the established 1e-6-Bohr derivative setting. All original physical gates remain:
stationary gradient, finite/symmetric Hessian, minimum and full finite-gradient
agreement. The replay additionally tests the two worst CDS-contact columns at
two displacement sizes. Its records remain explicitly isolated and ineligible
for production until independent review and source-bound admission. The original
failure is retained; neither acceptance thresholds nor scientific method changed.

Raw Hessian/CDS arrays stay on shared storage under this directory and
`calcs/dft/hyp_cata_adi_ethyl/neu/native_thermal/2d245a663438d6ee/`.
Their locations and byte hashes are in the frozen receipt and matrix review.

The full finer-step replay has now passed: projected minimum +14.53 cm⁻¹,
maximum asymmetry 7.11e-12 Ha/Bohr², and all four worst-contact finite-gradient
errors below 0.014%. An [independent review](finer_independent_review.json)
passes 223 checks, including the raw spectrum, stored gradients and agreement
with the independently reconstructed matrix. A [fully qualified cache](admission/preparation.json)
has been prepared with the original geometry, result metadata and active energy
bytes unchanged. Canonical admission is [waiting for the other molecule writers](admission/adoption.json);
the original failed state is still preserved. No canonical correction is claimed yet.

The cache has now been [adopted](admission/adoption.json) after exclusive writer
locks became available. Both readers accept the canonical correction; geometry,
raw result metadata and active-energy bytes are unchanged. The original failed
state is preserved in `admission/superseded_canonical/`, and the verified prepared
snapshot remains retained byte-for-byte. Repeating the admission command verifies
the existing adoption and archive without writing again.
