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
        family="pyridine-multi-e",
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
]
