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
        # MV2+/+. and MV+./0 are textbook; values vs Fc/Fc+ in MeCN — VERIFY vs OROP.
        events=[dict(event="ox2->ox1", exp_V_vs_Fc=-0.45, note="approx, verify"),
                dict(event="ox1->neu", exp_V_vs_Fc=-0.88, note="approx, verify")],
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
