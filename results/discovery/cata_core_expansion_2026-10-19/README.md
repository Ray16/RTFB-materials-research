# Alternative cata-annulated imide cores

The [source study](https://doi.org/10.3762/bjoc.20.155) reports four bis-phenyl/bis-hexyl
parents. Their eight CH2Cl2/Fc potentials remain source-parent measurements and
are not MeCN calibration labels or properties of either proposed graft.

[Source graphs](source_parents.json) are reconstructed from the deposited CIF bond
lists. The half-molecule 7-Ph requires its explicit crystallographic symmetry
mate. Two independent 8-Ph molecules give the same graph; their CH2Cl2 solvent
component is retained separately in the curation. The hexyl analogues follow the
source's declared N-cap replacement and NMR hydrogen counts. Four-carbonyl,
closed-shell, formula, symmetry and bond-distance checks pass. Generated
[parent drawings](parents.png) were compared with the source scheme.

| Prepared exact graft | Atoms including H | Matched control | Perturbation |
| --- | ---: | --- | --- |
| `hyp_cata_ndi_ethyl` | 48 | `design_ndi_ethyl` | Core-topology isomers with the same formula, N-ethyl cap and N-graft anchor |
| `hyp_cata_adi_ethyl` | 54 | `hyp_cata_ndi_ethyl` | Add one fused benzene; core size changes explicitly |

[Proposals](proposals.json), [drawings](prepared_grafts.png) and the
[isolated build receipt](embedding_receipt.json) preserve exact identities and
initial conformer ensembles. The shared candidate registry and manifest are
unchanged. Mixed N-ethyl/N-graft precursor isolation, selectivity and exact-graft
performance remain unverified. The declared two-electron ladders are hypotheses.

The [actual-start verification](uma_start_verification.json) separates resident
predictors and completed records from submitted controllers. The
[six-state UMA plan](uma_plan.json) and [launch record](uma_launch.json) track
gas preoptimization separately from source curation and solution properties.
The [completion review](uma_completion.json) verifies all six gas UMA states and
terminal workers, including memory headroom. No shared registration or solution
DFT calculation has occurred. The cata-NDI dianion has a 0.2348-eV UMA spin gap,
requiring a subsequent DFT alternate-spin comparison. Each task uses all initial conformers and the parity-consistent spin scan,
one reserved GPU, two assigned CPUs and capped native pools. One actual attempt
is permitted; failures preserve their outputs. Existing verified case records
can be inspected before any later resume. UMA energies never fill solution
potential or kinetic fields.

Raw CIFs and the scheme remain on shared storage under
`data/raw/discovery_sources/2026-10-19/cata_diimides/`; their publication URLs and
byte hashes are recorded in [sources.json](../../../data/raw/discovery_sources/2026-10-19/cata_diimides/sources.json).
Run `curate.py` with the recorded redox interpreter and native thread caps to
reconstruct the preparation. The [summary](summary.json) owns the check count.
The source catalog's 2026-10-18 cata entry describes the earlier retrieval stage;
this review completes its pending graph transcription.
