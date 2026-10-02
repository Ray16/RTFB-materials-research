"""Standalone (un-tethered) forms of the starting candidates, for the standalone-vs-grafted
2x2 (structure x environment). 'Standalone' = the redox core with the Merrifield benzyl handle
replaced by a minimal methyl cap; the two quinones keep their methoxy groups and are therefore
EXACTLY the Candidates.xlsx / D3TaLES molecules (enables a direct D3TaLES reorg cross-check).

Direct SMILES (no scaffold decoration). Same state schema as config/redox_groups.py:
states = (label, charge, mult, n_e); first state is the resting state.

NOTE: the standalone viologen (N,N'-dimethyl = methyl viologen) is already in the pipeline as
`methyl_viologen` (config/validation.py) and is NOT repeated here — the plot reuses it.
"""

STANDALONE = [
    dict(
        id="ethylviologen_sa",
        name="ethyl viologen (N,N'-diethyl-4,4'-bipyridinium)",
        # minimal capped analogue (the "methyl-capped" row of the 2x2; NOT the literal
        # pre-grafting precursor -- see the module docstring); NOT a screening
        # candidate. Tagged 'validation' to keep the two standalone viologens consistent.
        family="validation",
        smiles="CC[n+]1ccc(-c2cc[n+](CC)cc2)cc1",
        states=[("ox2", 2, 1, 0), ("ox1", 1, 2, -1), ("neu", 0, 1, -2)],
    ),
    dict(
        id="pmdi_sa",
        name="N,N'-dimethyl pyromellitic diimide",
        family="imide (n-type)",
        smiles="Cn1c(=O)c2cc3c(=O)n(C)c(=O)c3cc2c1=O",
        states=[("neu", 0, 1, 0), ("red1", -1, 2, -1), ("red2", -2, 1, -2)],
    ),
    dict(
        id="ndi_ammonium_sa",
        name="N-methyl-N'-(3-trimethylammoniopropyl) naphthalene diimide",
        family="imide (n-type)",
        smiles="CN1C(=O)c2ccc3c4c(ccc(c24)C1=O)C(=O)N(CCC[N+](C)(C)C)C3=O",
        states=[("ox", 1, 1, 0), ("red1", 0, 2, -1), ("red2", -1, 1, -2)],
    ),
    dict(
        id="mophquinone_sa",
        name="2-(4-methoxyphenyl)-1,4-benzoquinone (xlsx / D3TaLES)",
        family="quinone (n-type)",
        smiles="COc1ccc(C2=CC(=O)C=CC2=O)cc1",
        states=[("neu", 0, 1, 0), ("red1", -1, 2, -1), ("red2", -2, 1, -2)],
    ),
    dict(
        id="dmophquinone_sa",
        name="2-(2,5-dimethoxyphenyl)-5-methoxy-1,4-benzoquinone (xlsx / D3TaLES)",
        family="quinone (n-type)",
        smiles="COC1=CC(=O)C(c2cc(OC)ccc2OC)=CC1=O",
        states=[("neu", 0, 1, 0), ("red1", -1, 2, -1), ("red2", -2, 1, -2)],
    ),
    # --- batch 2: merrifield_multielectron_smiles.xlsx (see config/merrifield_multielectron.py)
    # Same convention: the 4-methylbenzyl resin tether replaced by a plain methyl cap.
    # The bis-tethered viologen's standalone analogue (both N capped with methyl) is
    # N,N'-dimethyl bipyridinium = `methyl_viologen`, already in config/validation.py —
    # NOT repeated here, the comparison reuses it.
    dict(
        id="aq_benzyloxy_sa",
        name="2-methoxyanthraquinone",
        family="quinone (n-type)",
        smiles="COc1ccc2c(c1)C(=O)c1ccccc1C2=O",
        states=[("neu", 0, 1, 0), ("red1", -1, 2, -1), ("red2", -2, 1, -2)],
    ),
    dict(
        id="nq_benzyloxy_sa",
        name="2-methoxy-1,4-naphthoquinone",
        family="quinone (n-type)",
        smiles="COC1=CC(=O)c2ccccc2C1=O",
        states=[("neu", 0, 1, 0), ("red1", -1, 2, -1), ("red2", -2, 1, -2)],
    ),
    dict(
        id="dtbc_phenol_sa",
        name="3,5-di-tert-butyl-2-methoxyphenol",
        family="phenol (p-type)",
        smiles="Oc1cc(C(C)(C)C)cc(C(C)(C)C)c1OC",
        states=[("neu", 0, 1, 0), ("ox", 1, 2, 1)],
    ),
    dict(
        id="aq_benzylamino_sa",
        name="2-(methylamino)anthraquinone",
        family="quinone (n-type)",
        smiles="CNc1ccc2c(c1)C(=O)c1ccccc1C2=O",
        states=[("neu", 0, 1, 0), ("red1", -1, 2, -1), ("red2", -2, 1, -2)],
    ),
]
