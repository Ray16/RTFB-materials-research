"""Second batch of Merrifield-grafted MULTI-ELECTRON candidates.

Source: data/raw/candidates/merrifield_multielectron_smiles.xlsx (13 rows, organized by
grafting chemistry: N-alkylation, Williamson ether, imidazole/Schiff-base metal linkers,
amine N-alkylation).

Same schema as redox_groups.py / starting_candidates.py: each core is a fragment carrying
one [*:1] dummy at the polymer-tether site; SCAFFOLD supplies the p-methylbenzyl resin
model. states = (label, charge, mult, n_e); first state is resting; n_e is electrons
transferred FROM the resting state (negative = reduction, positive = oxidation).

WHAT WAS TAKEN FROM THE SHEET AND WHAT WAS NOT
----------------------------------------------
Included (3 clean organics):
  - bisviologen    : row "4,4'-Bipyridine -> Viologen (xPS-V2+)". BOTH pyridine N's alkylate
                     resin CH2Cl, so the grafted model carries TWO benzyl tethers. This is
                     genuinely distinct from the mono-tethered `viologen` in
                     config/starting_candidates.py -- not a duplicate.
  - aq_benzyloxy   : row "2-Hydroxyanthraquinone", Williamson ether. Note the existing
                     `anthraquinone` entry in redox_groups.py is the CH2-O-CH2 linked
                     isomer; this is the direct aryl-O-CH2 ether, a different molecule.
  - nq_benzyloxy   : row "1,4-Naphthoquinol", grafted as lawsone (2-OH-1,4-NQ) O-benzyl
                     ether -- the sheet's own recommended route, no post-oxidation needed.

DEDUPED (already in the library, not rebuilt here):
  - "N-Monomethyl-4,4'-bipyridinium -> MV2+" is exactly `viologen` in starting_candidates.py
    (N-benzyl-N'-methyl bipyridinium). Skipped.

REPAIRED + FLAGGED (the sheet's chemistry does not close; see `flag` key):
  - dtbc_phenol    : the sheet grafts ONE catechol -OH of 3,5-di-tert-butylcatechol as a
                     benzyl ether and then claims the 2e catechol/o-quinone couple. That
                     couple needs BOTH oxygens -- a mono-O-benzyl ether is locked as an
                     ether and cannot form the o-quinone. What the grafted species can
                     actually do is a 1e phenol -> phenoxyl radical-cation oxidation, which
                     is what is modelled here. NOT a 2e centre as written.
  - aq_benzylamino : the sheet grafts 9-aminoanthracene through its C9 -NH2, then oxidizes
                     anthracene -> anthraquinone. But C9 is exactly where the quinone
                     carbonyl forms, so the oxidation destroys the tether. Modelled instead
                     as the 2-(benzylamino)anthraquinone regioisomer, where the amine sits
                     on a peripheral ring carbon and survives oxidation. Same 2e AQ couple,
                     C-N linkage as the sheet intends.

EXCLUDED (5 transition-metal rows: Fe/Mn-THPP porphyrin, Co/Mn-salen, Mo/W/Re-oxo Schiff
base). Three independent reasons, documented in FINDINGS.md:
  1. The sheet's "resin-attached" SMILES are placeholders, not the grafted complexes --
     Fe(THPP) and Mn(THPP) both give the H2THPP FREE BASE (no metal, identical string);
     WO2 gives the WO2(acac)2 precursor; Re(V) gives methyltrioxorhenium. Nothing to build.
  2. SMD is not parameterized for these metals. The MNSOL CDS term returns essentially the
     same value for Fe/Mn/Co/Mo/W (all within 0.002 eV on a fixed test geometry) because no
     element-specific surface tension exists for them -- and solvation dominates dG for
     charged metal complexes.
  3. The Mn/Mo/W/Re couples are catalytic oxo-transfer (M=O bond made/broken), so the
     Nelsen 4-point lambda_i -- which assumes one molecule on two adiabatic surfaces with no
     bond change -- is not a meaningful descriptor for them at any spin multiplicity.
  Fe(THPP) Fe3+/Fe2+ and Co(salen) Co2+/Co3+ ARE genuine outer-sphere 1e couples and could
  be run later as a flagged sub-study (needs metalated structures built by hand + spin-state
  scans); they are out of scope for this batch.
"""

SCAFFOLD = "Cc1ccc(C[*:1])cc1"

GROUPS = [
    dict(
        id="bisviologen",
        name="N,N'-bis(benzyl)-4,4'-bipyridinium (bis-tethered viologen)",
        family="pyridine-multi-e",
        frag="[*:1][n+]1ccc(-c2cc[n+](Cc3ccc(C)cc3)cc2)cc1",
        states=[
            ("ox2", +2, 1,  0),   # dication (resting)
            ("ox1", +1, 2, -1),   # 1e reduced -> radical cation
            ("neu",  0, 1, -2),   # 2e reduced -> neutral
        ],
    ),
    dict(
        id="aq_benzyloxy",
        name="2-(benzyloxy)anthraquinone",
        family="quinone (n-type)",
        frag="[*:1]Oc1ccc2c(c1)C(=O)c1ccccc1C2=O",
        states=[
            ("neu",  0, 1,  0),   # neutral quinone (resting)
            ("red1",-1, 2, -1),   # semiquinone radical anion
            ("red2",-2, 1, -2),   # dianion
        ],
    ),
    dict(
        id="nq_benzyloxy",
        name="2-(benzyloxy)-1,4-naphthoquinone (from lawsone)",
        family="quinone (n-type)",
        frag="[*:1]OC1=CC(=O)c2ccccc2C1=O",
        states=[
            ("neu",  0, 1,  0),
            ("red1",-1, 2, -1),
            ("red2",-2, 1, -2),
        ],
    ),
    dict(
        id="dtbc_phenol",
        name="3,5-di-tert-butyl-2-(benzyloxy)phenol (mono-O-benzyl DTBC)",
        family="phenol (p-type)",
        frag="Oc1cc(C(C)(C)C)cc(C(C)(C)C)c1O[*:1]",
        flag="1e-only: mono-O-benzyl ether cannot form the o-quinone; sheet's 2e "
             "catechol/o-quinone couple is not accessible on the grafted species. "
             "ALSO mechanistically out of scope: ArOH-+ is a strong acid (pKa <~ -2 in "
             "MeCN), so real phenol oxidation is EC (electron transfer then O-H "
             "deprotonation), not a reversible outer-sphere couple.",
        # Report the numbers, but do NOT rank it: the Nelsen 4-point lambda_i assumes one
        # molecule on two adiabatic surfaces with NO bond made or broken. O-H cleavage
        # violates that exactly as M=O formation does for the excluded metal-oxo rows, so
        # ranking this against genuine outer-sphere couples would not be like-for-like.
        rankable=False,
        states=[
            ("neu",  0, 1,  0),   # neutral phenol (resting)
            ("ox",  +1, 2, +1),   # 1e oxidation -> phenoxyl radical cation
        ],
    ),
    dict(
        id="aq_benzylamino",
        name="2-(benzylamino)anthraquinone (repaired 9-aminoanthracene route)",
        family="quinone (n-type)",
        frag="[*:1]Nc1ccc2c(c1)C(=O)c1ccccc1C2=O",
        flag="regiochemistry repaired: sheet's C9-NH2 tether is destroyed by the "
             "anthracene->anthraquinone oxidation (C9 becomes the carbonyl); modelled at C2",
        states=[
            ("neu",  0, 1,  0),
            ("red1",-1, 2, -1),
            ("red2",-2, 1, -2),
        ],
    ),
]
