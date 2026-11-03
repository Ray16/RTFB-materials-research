# Exact 2-cyanophenothiazine precursor access

One new direct N-linked design now has reviewed exact precursor-listing access:
`proposal_ptz_176c4d2f397734ea`. The [independent check](independent_check_03.json)
matches its N-H precursor to TCI P2308 (CAS 38642-74-9), using the manufacturer's
[primary listing](https://www.tcichemicals.com/US/en/p/P2308), deposited PubChem
SID 468592587 and CID 98841. Its canonical graph and InChIKey agree. Other cyano
positions, fluorinated precursors and N-hydroxyethyl/ether-linked designs do not
inherit this evidence. A listing is not acquired material, lot QC or an
independent synthesis; stock and price are not asserted here.

Cao et al., *Crystal Growth & Design* 2024, DOI
[10.1021/acs.cgd.4c00789](https://doi.org/10.1021/acs.cgd.4c00789), report the
N-butyl molecular analogue in their official Supporting Information, p.2.
The public ACS/Figshare PDF was consulted directly, including its reported
NMR/HRMS assignments. NaH/DMF N-alkylation preserves the cyano/phenothiazine
connectivity in that molecule. This is a molecular analogue, not a demonstration
of benzyl attachment to crosslinked beads or retained battery redox function.
The raw PDF, primary metadata and retrieval records remain on shared storage;
[source_review.json](source_review.json) binds their paths and hashes. Only small
identity/DOI records and our interpretation are versioned here.

The [updated exact evidence](reviewed_candidate_evidence.json) changes precursor
access only. The [prospective gate](prospective_gate.json) still holds the full
redox campaign for exact-graft selectivity, a balanced actual-condition model
and redox retention. Proposed NaH/DMF bead conditions remain explicitly
transferred from the separate unsubstituted bead report; no actual-condition G
or yield is calculated. No operating property transfers from either source.

The maintained merger published a fresh complete pool at
`results/candidates/pool/2026-11-02_ptz_precursor_access/`; the stable
`results/candidates/pool/current.json` selects it. All 912 graphs, prior aliases,
declarations and computed evaluations are preserved. There are still three
route passes, 909 holds and zero complete material admissions. The earlier
[direct-phenothiazine route/pool review](../phenothiazine_postfunctionalization_2026-11-02/README.md)
and frozen manuscript epoch remain intact. Four identity tests reject wrong
supplier/CID/structure identifiers; they are not quantum or experimental tests.

Recheck with `python scripts/validation/check_ptz_precursor_access.py --help`.
The independent check also verifies that curated metadata are byte-identical
to the consulted primary records and that no previous property changed. The
first independent check and its exact review-source bytes are preserved under
`review_sources/`; the second adds curated-copy bindings, and the third checks the exact previous catalog as well.
