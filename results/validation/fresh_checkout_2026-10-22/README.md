# Clean-checkout CPU reproduction

The [checked record](summary.json) verifies the small example from a new sparse
Git checkout of commit `c09bc544a7946751fb31a2f429c4d83c80344401`, using a newly
created Python 3.11 environment. Both system and user site packages were disabled.
The project loaded from that checkout, and every scientific dependency loaded
from the new environment. The dependency check found no broken requirements.

The example reproduced 32 recipes / 17 exact designs, matched every descriptor
with serial replay, held all seventeen designs for missing route evidence and
submitted no quantum jobs. The regenerated conditional-grafting CSV matches
exactly; its rendered RGB pixels also match the checked figure exactly.

Inspect [the example output](smoke.json), [environment isolation](isolation.json),
[package/download identities](installation.json), [full versions](pip_freeze.txt)
and [dependency check](pip_check.txt). Raw checkout, installation logs, complete
wheel report and regenerated outputs remain at
`calcs/release_validation/2026-10-22_fresh_checkout/`.

This is a real clean-source and fresh-environment CPU reproduction. It uses
curated scientific data to regenerate a figure; it does not rerun quantum
energies. The checkout contains no raw DFT/UMA calculation trees. Local shared
Git objects avoid duplicating history but do not constitute a public archival
release. Full quantum/CUDA reproduction and quantum interruption equivalence
remain separate requirements.
