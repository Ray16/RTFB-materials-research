# DIPEA basicity: source origin precedes validation

[summary.json](summary.json) records the checked source notes, exact DIPEA identity,
component hashes and model sensitivities. [source_spec.json](source_spec.json)
preserves the source specification. The pipeline entry point is
[review_basicity_origins.py](../../../scripts/pipeline/review_basicity_origins.py);
the reusable source policy is [basicity_reference_origins.json](../../../config/basicity_reference_origins.json).

The author-posted [Heldebrant et al. source](https://www.pnnl.gov/sites/default/files/media/file/at/publications/322_427201153539.pdf)
lists Hünig's base at pKaH 18.6 in MeCN with a specific PCM-estimate note in
Table 1, while Table 2 uses a generic experimental label. We preserve this
origin conflict and exclude the value from independent measured validation.
The subsequent Mannisto et al. citation points back to that source; an
independent experimental determination was not established by this audit.

At the stated comparison temperature, using that computed literature value
would give a capture term of −47.38 kJ/mol. With the stored HCl-reference sums,
the two BQ sensitivities become +10.26/+8.66 kJ/mol, compared with the all-DFT
DIPEA values −13.34/−14.94. This is a model-to-model sensitivity, not a
calibrated prediction or uncertainty interval. HCl cancels in the complete
all-DFT reaction, so replacing just its capture term changes shared-error
cancellation. No primary reaction energy or candidate decision was changed.

Third-party PDFs and rendered pages remain on shared storage under
`calcs/literature/2026-10-23_dipea_reference_origin`; their hashes are recorded.
They are not redistributed in this code checkpoint. The independently measured
DIPEA MeCN anchor, source activity/temperature conventions and complete-route
validation remain unresolved.
