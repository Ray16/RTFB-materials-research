"""Starting candidates — a focused set of promising MULTI-CHARGE (multi-electron) redox
actives, grafted onto the Merrifield benzylic site exactly like config/redox_groups.py.

Source: starting_candidates/Candidates.xlsx (literature non-aqueous flow-battery actives +
two D3TaLES methoxy-quinone exemplars). Each core is expressed as a fragment carrying one
[*:1] dummy at the polymer-tether site (the benzyl handle from SCAFFOLD attaches there);
the other substituent site on the symmetric diimide/bipyridinium cores is capped with the
paper's group (methyl / trimethylammoniopropyl) so the redox core is faithfully represented.

Same schema as redox_groups.py: states = list of (label, charge, mult, n_e); first state is
the resting state; n_e = electrons transferred from resting (negative = reduction).

Notes / provenance:
- methylviologen from the sheet == the existing `viologen` entry (N-benzyl-N'-methyl
  bipyridinium) in redox_groups.py — DEDUPED, not repeated here.
- ethylviologen: N'-ethyl analogue of that viologen (the sheet's ethylviologen).
- pmdi: pyromellitic diimide anolyte (sheet: N-2-pentyl PMDI); the floppy 2-pentyl (an
  unassigned stereocenter) and PEG chain are polymer-model artifacts, so the tether N is
  the graft and the other N is methyl-capped.
- ndi_ammonium: ammonium-functionalized naphthalene diimide; one imide N is the benzyl
  tether, the other keeps the paper's charge-stabilizing 3-(trimethylammonio)propyl group
  (permanent +1 on the resting core -> net +1/0/-1 across the 2e reduction).
- mophquinone / dmophquinone: the two D3TaLES methoxy-benzoquinones; one methoxy is
  promoted to a benzyloxy tether (alkoxy->alkoxy, ~isoelectronic) to graft the core.
"""

SCAFFOLD = "Cc1ccc(C[*:1])cc1"

GROUPS = [
    dict(
        id="ethylviologen",
        name="N-benzyl-N'-ethyl-4,4'-bipyridinium (ethyl viologen)",
        family="pyridine-multi-e",
        frag="[*:1][n+]1ccc(-c2cc[n+](CC)cc2)cc1",
        states=[
            ("ox2", +2, 1,  0),   # dication (resting)
            ("ox1", +1, 2, -1),   # 1e reduced -> radical cation
            ("neu",  0, 1, -2),   # 2e reduced -> neutral
        ],
    ),
    dict(
        id="pmdi",
        name="N-benzyl-N'-methyl pyromellitic diimide",
        family="imide (n-type)",
        frag="[*:1]n1c(=O)c2cc3c(=O)n(C)c(=O)c3cc2c1=O",
        states=[
            ("neu",  0, 1,  0),   # neutral (resting)
            ("red1",-1, 2, -1),   # radical anion
            ("red2",-2, 1, -2),   # dianion
        ],
    ),
    dict(
        id="ndi_ammonium",
        name="N-benzyl N'-(3-trimethylammoniopropyl) naphthalene diimide",
        family="imide (n-type)",
        frag="[*:1]N1C(=O)c2ccc3c4c(ccc(c24)C1=O)C(=O)N(CCC[N+](C)(C)C)C3=O",
        states=[
            ("ox",  +1, 1,  0),   # resting: NDI core neutral, +1 from the ammonium tether
            ("red1", 0, 2, -1),   # NDI radical anion (net neutral)
            ("red2",-1, 1, -2),   # NDI dianion (net -1)
        ],
    ),
    dict(
        id="mophquinone",
        name="2-(4-benzyloxyphenyl)-1,4-benzoquinone",
        family="quinone (n-type)",
        frag="[*:1]Oc1ccc(C2=CC(=O)C=CC2=O)cc1",
        states=[
            ("neu",  0, 1,  0),   # neutral quinone (resting)
            ("red1",-1, 2, -1),   # semiquinone radical anion
            ("red2",-2, 1, -2),   # hydroquinone dianion
        ],
    ),
    dict(
        id="dmophquinone",
        name="2-(2-methoxy-5-benzyloxyphenyl)-5-methoxy-1,4-benzoquinone",
        family="quinone (n-type)",
        frag="COC1=CC(=O)C(c2cc(O[*:1])ccc2OC)=CC1=O",
        states=[
            ("neu",  0, 1,  0),
            ("red1",-1, 2, -1),
            ("red2",-2, 1, -2),
        ],
    ),
]
