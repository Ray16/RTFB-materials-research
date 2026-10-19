# Solution-relaxation evidence and cache contracts

Implementation controls are now separate from per-case verification. A new summary
defaults to unverified; producers may assert a passed implementation control only
after checking bound CPU/independent-hcore and requested-backend evidence. That
control does not independently replay every target calculation or validate a rate.
The paired runner verifies supplied control hashes before launching SCF work.

Both documented SCF cache schemas are recognized. Reuse requires matching geometry,
charge, spin, solvent response, frozen reference density and numerical source hashes.
A density must have its own byte binding or an explicit source-bound independent
review of both JSON and NPZ. Shapes, symmetry, finite values, electron traces,
energy units and open-shell spin diagnostics are checked. Unverified pool entries
are preserved and skipped; inconsistent records in the current output are refused.
Fresh records bind both density artifact and matrix hashes.

The [archived pilot review](pilot_cache_review.json) rechecks all six preserved
PTZ densities and independently recomputes eight directional energy differences
and same-geometry CDS cancellation: **16 checks, zero failures**. The symmetric
mean remains 1.126525 eV. No quantum SCF was repeated and no raw file was changed.
The original producer source was recovered by exact hash from Git and retained
under `original_producers/`, preserving historical evidence after the code fix.

[Qualification metadata](qualification.json) separately tracks endpoint-model
and symmetric-model eligibility for four existing controls. These results remain
composite electronic solution descriptors, not experimentally validated rates.
The flexible PTZ endpoints require the reciprocal local-basin review in
`results/validation/ptz_basin_association_2026-10-19/` before its conformer flag can
be cleared. Active-SP stationarity and a complete conformer ensemble remain separate.
