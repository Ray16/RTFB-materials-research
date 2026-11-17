"""Redox-group library for decorating the Merrifield monomer's benzylic site.

The Cl of the chloromethyl-polystyrene monomer is replaced by each group below.
The scaffold is a truncated molecular model of the resin repeat unit: a para-alkyl
benzyl handle (the alkyl stub stands in for the polymer backbone), with a dummy atom
[*:1] marking the benzylic carbon where the redox group attaches (i.e. where Cl was).

Each group defines the redox STATES relevant to it, with explicit total `charge` and
spin `mult` (multiplicity) — UMA needs both. `n_e` on a state is electrons transferred
from the resting state (negative = reduction, positive = oxidation).
"""

# Truncated resin model: 4-methylbenzyl handle. The para-CH3 is the minimal backbone
# stub (an alkyl donor para to the redox site); its small electronic effect cancels in
# group-to-group comparisons. For a fidelity check, swap CH3 -> CH(C)C (isopropyl),
# which reproduces the secondary backbone methine. [*:1] marks the benzylic C (ex-Cl).
SCAFFOLD = "Cc1ccc(C[*:1])cc1"

# Each group: attachment fragment carries one [*:1] dummy bonded to the attaching atom.
# states: list of (label, charge, mult, n_e). The first state is the resting state.
GROUPS = [
    dict(
        id="pyridinium",
        name="N-benzylpyridinium",
        family="pyridine",
        frag="[*:1][n+]1ccccc1",
        states=[
            ("ox",  +1, 1,  0),   # resting cation
            ("red",  0, 2, -1),   # 1e reduction -> neutral radical
        ],
    ),
    dict(
        id="cyanopyridinium",
        name="N-benzyl-4-cyanopyridinium",
        family="pyridine",
        frag="[*:1][n+]1ccc(C#N)cc1",
        states=[
            ("ox",  +1, 1,  0),
            ("red",  0, 2, -1),
        ],
    ),
    dict(
        id="viologen",
        name="N-benzyl-N'-methyl-4,4'-bipyridinium (viologen)",
        family="pyridine-multi-e",
        frag="[*:1][n+]1ccc(-c2cc[n+](C)cc2)cc1",
        states=[
            ("ox2", +2, 1,  0),   # dication (resting)
            ("ox1", +1, 2, -1),   # 1e reduced -> radical cation
            ("neu",  0, 1, -2),   # 2e reduced -> neutral (closed shell)
        ],
    ),
    dict(
        id="phenothiazine",
        name="N-benzylphenothiazine",
        family="amine (p-type)",
        frag="[*:1]N1c2ccccc2Sc2ccccc21",
        states=[
            ("neu",  0, 1,  0),   # neutral (resting)
            ("ox",  +1, 2, +1),   # 1e oxidation -> radical cation
        ],
    ),
    dict(
        id="tempo",
        name="4-(benzyloxy)-TEMPO",
        family="nitroxide",
        frag="[*:1]OC1CC(C)(C)N([O])C(C)(C)C1",
        states=[
            ("rad",  0, 2,  0),   # nitroxide radical (resting)
            ("ox",  +1, 1, +1),   # -> oxoammonium
            ("red", -1, 1, -1),   # -> hydroxylamine anion
        ],
    ),
    dict(
        id="anthraquinone",
        name="2-(benzyloxymethyl)anthraquinone",
        family="quinone (n-type)",
        frag="[*:1]OCc1ccc2C(=O)c3ccccc3C(=O)c2c1",
        states=[
            ("neu",  0, 1,  0),   # neutral (resting)
            ("red1",-1, 2, -1),   # radical anion
            ("red2",-2, 1, -2),   # dianion
        ],
    ),
]
