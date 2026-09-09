"""Validation set: parent redox cores with KNOWN experimental redox potentials in MeCN,
used to test the pipeline against measurement before trusting decorated-monomer rankings
(docs/PLAN.md §V, hard gate). Same schema as config/redox_groups.py so they run through
the identical pipeline (build -> UMA -> DFT+SMD -> redox).

`exp_V_vs_Fc` on an event is the experimental potential (V vs Fc/Fc+ in MeCN). These are
literature ANCHORS with citations and are marked provisional — the primary, consistently
referenced experimental set is OROP/ROP313 (data/raw/validation/SI_data_redox_paper),
matched to these molecules by structure. Verify/replace anchors against OROP + ReSolvedDB
before using them to calibrate. Do NOT fit away discrepancies — diagnose the physics.

Ferrocene is the internal reference (E° vs Fc/Fc+ = 0 by definition); it also needs a
metallocene geometry that RDKit cannot embed, so it is flagged special (build separately).

REFERENCE-FRAME AUDIT (2026-09): each anchor was checked against the established MeCN-vs-Fc
window for its couple. TEMPO+/TEMPO (+0.24), phenothiazine+./PTZ (+0.26), 9,10-anthraquinone
(-1.28/-1.90), and N-methylpyridinium (-1.8) all sit inside their vs-Fc windows -> correctly
labeled. ONLY methyl viologen was wrong: -0.45/-0.88 are MeCN vs SCE, mislabeled vs Fc, now
converted (see below). Remaining TODO: back every anchor with a primary MeCN CV citation.
"""

VALIDATION = [
    dict(
        id="ferrocene",
        name="ferrocene (internal reference)",
        smiles="[Fe].c1ccc[cH-]1.c1ccc[cH-]1",   # metallocene: needs special geometry
        special_geometry=True,
        states=[("neu", 0, 1, 0), ("ox", 1, 2, +1)],
        events=[dict(event="ox->neu", exp_V_vs_Fc=0.00, note="defines the Fc/Fc+ scale")],
    ),
    dict(
        id="methyl_viologen",
        name="methyl viologen (N,N'-dimethyl-4,4'-bipyridinium)",
        smiles="C[n+]1ccc(-c2cc[n+](C)cc2)cc1",
        states=[("ox2", 2, 1, 0), ("ox1", 1, 2, -1), ("neu", 0, 1, -2)],
        # CORRECTED 2026-09 (reference-frame audit): the prior -0.45/-0.88 were MeCN vs SCE
        # (docs/PLAN.md states them "vs SCE"), MISLABELED here as vs Fc. They sit ~0.4 V above
        # the established MeCN-vs-Fc viologen window (1st red ~-0.8, 2nd red ~-1.25), i.e. in the
        # vs-SCE window. Converted to Fc via Pavlishchuk & Addison (Inorg. Chim. Acta 2000):
        # Fc/Fc+ = +0.40 V vs SCE in MeCN  =>  E_vs_Fc = E_vs_SCE - 0.40.
        # PROVISIONAL (+/-~0.05 V); replace with a primary MeCN CV vs Fc when available. The
        # reference-free wave spacing (E1-E2 ~0.44 V) independently matches experiment and our
        # computed ladder, so the shift is a referencing fix, not a fit to our numbers.
        events=[dict(event="ox2->ox1", exp_V_vs_Fc=-0.85,
                     note="MeCN vs SCE (-0.45) - 0.40 (Pavlishchuk-Addison); provisional"),
                dict(event="ox1->neu", exp_V_vs_Fc=-1.28,
                     note="MeCN vs SCE (-0.88) - 0.40 (Pavlishchuk-Addison); provisional")],
    ),
    # --- Explicit PF6- ion-pair species for the viologen fix (released-counterion scheme,
    # docs/PLAN.md). Each keeps its NATURAL number of PF6- so every assembly is NEUTRAL
    # (best for continuum SMD); each reduction releases one free PF6-. E° is assembled from
    # these by src/redox/ionpair.py as a cross-species reaction, NOT redox.py's adjacent-
    # charge pairing — so each species carries a SINGLE state (no auto-couples). The SMILES
    # formal charges only seed RDKit; the per-state (charge, mult) sets the actual DFT charge
    # (UMA re-scans spin). PF6- = F[P-](F)(F)(F)(F)F.
    dict(
        id="pf6",
        name="hexafluorophosphate (free PF6- anion)",
        smiles="F[P-](F)(F)(F)(F)F",
        states=[("anion", -1, 1, 0)],   # 70 e- (even) -> singlet
    ),
    dict(
        id="mv_ip2",
        name="methyl viologen dication . 2 PF6- (ion pair, net 0)",
        smiles="C[n+]1ccc(-c2cc[n+](C)cc2)cc1.F[P-](F)(F)(F)(F)F.F[P-](F)(F)(F)(F)F",
        states=[("s0", 0, 1, 0)],       # MV2+ + 2 PF6- ; net 0, 238 e- (even) -> singlet
    ),
    dict(
        id="mv_ip1",
        name="methyl viologen radical-cation . 1 PF6- (ion pair, net 0)",
        smiles="C[n+]1ccc(-c2cc[n+](C)cc2)cc1.F[P-](F)(F)(F)(F)F",
        states=[("s0", 0, 2, 0)],       # MV+. + 1 PF6- ; net 0, 169 e- (odd) -> doublet
    ),
    dict(
        id="tempo_parent",
        name="TEMPO (2,2,6,6-tetramethylpiperidine-1-oxyl)",
        smiles="CC1(C)CCCC(C)(C)N1[O]",
        states=[("rad", 0, 2, 0), ("ox", 1, 1, +1), ("red", -1, 1, -1)],
        events=[dict(event="ox->rad", exp_V_vs_Fc=+0.24, note="TEMPO+/TEMPO, approx verify")],
    ),
    dict(
        id="phenothiazine_parent",
        name="10H-phenothiazine",
        smiles="c1ccc2c(c1)Nc1ccccc1S2",
        states=[("neu", 0, 1, 0), ("ox", 1, 2, +1)],
        events=[dict(event="ox->neu", exp_V_vs_Fc=+0.26, note="PTZ+./PTZ, approx verify")],
    ),
    dict(
        id="anthraquinone_parent",
        name="9,10-anthraquinone",
        smiles="O=C1c2ccccc2C(=O)c2ccccc21",
        states=[("neu", 0, 1, 0), ("red1", -1, 2, -1), ("red2", -2, 1, -2)],
        events=[dict(event="neu->red1", exp_V_vs_Fc=-1.28, note="AQ/AQ-. approx verify"),
                dict(event="red1->red2", exp_V_vs_Fc=-1.90, note="approx verify")],
    ),
    dict(
        id="methylpyridinium",
        name="N-methylpyridinium",
        smiles="C[n+]1ccccc1",
        states=[("ox", 1, 1, 0), ("red", 0, 2, -1)],
        events=[dict(event="ox->red", exp_V_vs_Fc=-1.8, note="hard to reduce, approx verify")],
    ),
]
