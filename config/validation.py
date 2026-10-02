"""Validation set: parent redox cores with KNOWN experimental redox potentials in MeCN,
used to test the pipeline against measurement before trusting decorated-monomer rankings
(docs/PLAN.md §V, hard gate). Same schema as config/redox_groups.py so they run through
the identical pipeline (build -> UMA -> DFT+SMD -> redox).

`exp_V_vs_Fc` on an event is the experimental potential (V vs Fc/Fc+ in MeCN). These are
literature ANCHORS with citations and are marked provisional — the primary, consistently
referenced experimental set is the ORGANIC OROP subset (193 systems, H/C/N/O/S/halogens,
MeCN+DMF, charges -2..+1). NOMENCLATURE (Neugebauer/Liu ROP313 paper): ROP313 = OROP (193
organic) + OMROP (120 ORGANOMETALLIC, transition-metal complexes, charges -4..+3). Our
pipeline validates against OROP ONLY (systems 1..193); the OMROP metal complexes are out
of scope. The clone at data/raw/validation/OROP/ holds the full ROP313 CSV (313 rows) but
only rows 1..193 are organic. Verify/replace these hand anchors against OROP experimental
values; do NOT fit away discrepancies — diagnose the physics.

Ferrocene is the internal reference (E° vs Fc/Fc+ = 0 by definition); it also needs a
metallocene geometry that RDKit cannot embed, so it is flagged special (build separately).

GROUNDING AUDIT (2026-10-02): every event now carries `grounded`. True only when the value
was read in the primary source (quoted in the note); False otherwise. Only grounded events
enter accuracy statistics / sigma derivation. Status: ferrocene (definition) and TEMPO
(Gerken & Stahl 2015) grounded; methyl viologen absolute values depend on an unverified
Ag/Ag+ -> Fc constant (its reference-free wave spacing is grounded, exp_dE12_V);
10H-phenothiazine (+0.26), 9,10-anthraquinone (-1.28/-1.90) and N-methylpyridinium (-1.8)
have NO verified primary citation. The sourced MeCN E1/E2 benchmark that supersedes these
for quinones/imides/viologens is config/benchmark_mecn.py
(data/raw/validation/two_wave_mecn/).
"""

VALIDATION = [
    dict(
        id="ferrocene",
        name="ferrocene (internal reference)",
        smiles="[Fe].c1ccc[cH-]1.c1ccc[cH-]1",   # metallocene: needs special geometry
        special_geometry=True,
        states=[("neu", 0, 1, 0), ("ox", 1, 2, +1)],
        events=[dict(event="ox->neu", exp_V_vs_Fc=0.00,
                     grounded=True,
                     note="Defines the scale: E(Fc+/Fc) = 0 vs Fc/Fc+ by construction.")],
    ),
    dict(
        id="methyl_viologen",
        family="pyridine-multi-e",
        name="methyl viologen (N,N'-dimethyl-4,4'-bipyridinium)",
        smiles="C[n+]1ccc(-c2cc[n+](C)cc2)cc1",
        states=[("ox2", 2, 1, 0), ("ox1", 1, 2, -1), ("neu", 0, 1, -2)],
        # GROUNDING AUDIT 2026-10-02. The previous -0.85/-1.28 V came from Bird & Kuhn's AQUEOUS
        # MV2+/+. couple (-0.446 V vs NHE) treated as "~ -0.45 vs SCE" in MeCN — a wrong
        # conversion that matched MeCN data only by coincidence. Primary MeCN source now:
        # Cook et al., ChemElectroChem 2017 (doi 10.1002/celc.201600536, open copy PMC5467523),
        # Table 1: E_a = -759 mV, E_b = -1179 mV vs Ag/Ag+ (10 mM, MeCN), 0.1 M TBAPF6, Pt, 0.1 V/s.
        # The paper reports no Fc calibration, so the vs-Fc values need a literature Ag/Ag+ -> Fc
        # constant we could not verify (-0.082 V, Pavlishchuk & Addison 2000 as transcribed
        # secondarily) -> grounded=False for the absolute values. The text states "the second
        # reduction is not as well behaved" -> E2 is not used. The wave spacing
        # E_a - E_b = 0.420 V is reference-free and grounded (exp_dE12_V).
        # N-benzyl-N'-methyl viologen (our grafted core) is measured DIRECTLY vs Fc at
        # -0.785 / -1.204 V by Kim et al., Chem. Eur. J. 2022 (doi 10.1002/chem.202200149,
        # PMC9310624; previously misattributed here to Cook 2017) — in config/benchmark_mecn.py.
        exp_dE12_V=0.420, exp_dE12_note="Cook 2017 Table 1: -759 - (-1179) mV, reference-free",
        events=[dict(event="ox2->ox1", exp_V_vs_Fc=-0.841, grounded=False,
                     note="Cook 2017 Table 1 -0.759 V vs Ag/Ag+(10 mM) - 0.082 (UNVERIFIED "
                          "secondary Ag/Ag+ -> Fc constant)"),
                dict(event="ox1->neu", exp_V_vs_Fc=-1.261, grounded=False,
                     note="Cook 2017 Table 1 -1.179 V vs Ag/Ag+ - 0.082 (unverified constant); "
                          "source: second reduction 'not as well behaved' -> not used")],
    ),
    # --- Explicit PF6- ion-pair species for the viologen fix (released-counterion scheme,
    # docs/PLAN.md). Each keeps its NATURAL number of PF6- so every assembly is NEUTRAL
    # (best for continuum SMD); each reduction releases one free PF6-. E° is assembled from
    # these by archive/src/redox/ionpair.py as a cross-species reaction, NOT redox.properties.potentials' adjacent-
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
        family="nitroxide",
        name="TEMPO (2,2,6,6-tetramethylpiperidine-1-oxyl)",
        smiles="CC1(C)CCCC(C)(C)N1[O]",
        states=[("rad", 0, 2, 0), ("ox", 1, 1, +1), ("red", -1, 1, -1)],
        events=[dict(event="ox->rad", exp_V_vs_Fc=+0.249, grounded=True,
                     note="Gerken & Stahl, ACS Cent. Sci. 2015 (doi 10.1021/acscentsci.5b00163, "
                          "PMC4827547): 'reversible nitroxyl/oxoammonium redox process at "
                          "E1/2 = 249 mV vs Fc/Fc+' in CH3CN (verified in text). Electrolyte "
                          "0.5 M KPF6, not a tetraalkylammonium salt.")],
    ),
    dict(
        id="phenothiazine_parent",
        family="amine (p-type)",
        name="10H-phenothiazine",
        smiles="c1ccc2c(c1)Nc1ccccc1S2",
        states=[("neu", 0, 1, 0), ("ox", 1, 2, +1)],
        events=[dict(event="ox->neu", exp_V_vs_Fc=+0.26, grounded=False,
                     note="REF (verify): Connelly & Geiger, Chem. Rev. 1996, 96, 877 / Fu et al. "
                          "JACS 2005, 127, 7227. SHAKY: N-H parent radical cation DEPROTONATES "
                          "(chemically irreversible) -> +0.26 better describes an N-ALKYL "
                          "phenothiazine; our candidate is N-benzyl (reversible), so anchor on "
                          "N-methylphenothiazine instead.")],
    ),
    dict(
        id="anthraquinone_parent",
        family="quinone (n-type)",
        name="9,10-anthraquinone",
        smiles="O=C1c2ccccc2C(=O)c2ccccc21",
        states=[("neu", 0, 1, 0), ("red1", -1, 2, -1), ("red2", -2, 1, -2)],
        events=[dict(event="neu->red1", exp_V_vs_Fc=-1.28, grounded=False,
                     note="REF (verify exact): Fu, Liu, Yu, Wang, Guo, JACS 2005, 127, 7227 "
                          "(270 organics, MeCN redox potentials) or Connelly & Geiger, Chem. "
                          "Rev. 1996, 96, 877. Standard aprotic AQ/AQ-. value; ~0.1 V uncertain."),
                dict(event="red1->red2", exp_V_vs_Fc=-1.90, grounded=False,
                     note="REF (verify exact): Fu et al. JACS 2005, 127, 7227. AQ-./AQ2-; ~0.1 V.")],
    ),
    dict(
        id="methylpyridinium",
        family="pyridine",
        name="N-methylpyridinium",
        smiles="C[n+]1ccccc1",
        states=[("ox", 1, 1, 0), ("red", 0, 2, -1)],
        events=[dict(event="ox->red", exp_V_vs_Fc=-1.8, grounded=False,
                     note="WEAKEST, UNVERIFIED: reduction is IRREVERSIBLE near the MeCN cathodic "
                          "limit -> -1.8 is an ESTIMATE (cf. Kosower pyridinyl radicals, JACS; "
                          "or Fu et al. JACS 2005, 127, 7227). NOT among our candidate families "
                          "(pyridine single-e) -> RECOMMEND drop or de-weight.")],
    ),
]
