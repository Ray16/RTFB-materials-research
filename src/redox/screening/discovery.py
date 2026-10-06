"""Systematic identification of graftable multi-electron redox candidates.

Turns any list of molecules (a database dump, a literature list, an enumeration) into the
same grafted-candidate rows the pipeline screens, with every rule explicit:

  1. FAMILY        substructure patterns (FAMILIES). Each family carries whether it is a
                   genuine 2e- (two resolved 1e- waves) motif and the evidence for that:
                   measured two-wave spacings in data/raw/validation/two_wave_mecn/ (tier A,
                   config/benchmark_mecn.py) or the source's own statement.
  2. GRAFTING      the Merrifield (chloromethyl-polystyrene) attachment chemistry as RDKit
                   reactions (GRAFTS): every nucleophilic site that can displace the benzylic
                   chloride gives one grafted product, modelled with the same 4-methylbenzyl
                   tether proxy as every current candidate (config/starting_candidates.py
                   SCAFFOLD "Cc1ccc(C[*:1])cc1").
  3. RULES         fully specified stereochemistry (CLAUDE.md), allowed elements, single
                   fragment.
  4. DESCRIPTORS   computed IDENTICALLY for pool and current candidates: monomer MW, nominal
                   specific capacity n F / (3.6 MW) with n = family electrons (nominal: the
                   pipeline's capacity uses the computed accessible path), Ertl SA of the
                   grafted species (the scorer input prep of redox.properties.
                   capacity_and_proxies).

`rediscover()` is the correctness test: run on the PRECURSORS of the current candidates,
the rules must regenerate exactly the grafted candidates the pipeline already screens.
"""
from __future__ import annotations

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, Descriptors

RDLogger.DisableLog("rdApp.*")

TETHER = "Cc1ccc(C[*:1])cc1"     # 4-methylbenzyl proxy of the polymer tether

# name -> (SMARTS, n_electrons, multi_e, evidence)
FAMILIES = {
    "p-quinone": ("[#6;R]1(=O)[#6;R]=,:[#6;R][#6;R](=O)[#6;R]=,:[#6;R]1", 2, True,
                  "two reversible waves in MeCN, spacing 0.79-0.94 V (benzoquinones, tier A)"),
    "o-quinone": ("[#6;R]1(=O)[#6;R](=O)[#6;R]=,:[#6;R][#6;R]=,:[#6;R]1", 2, True,
                  "two quasi-reversible waves, spacing 0.41 V (1,2-naphthoquinone, tier A)"),
    "bis-imide": (None, 2, True,
                  "two waves, spacing 0.61-0.64 V (PMDIs, tier A) / 0.44 V (NDI)"),
    "mono-imide": ("O=[#6;R](c)[#7;R][#6;R](=O)c", 1, False,
                   "single reversible wave only (N-alkyl phthalimides, Daub 2021)"),
    "viologen": ("[n+;R]1ccc(-c2cc[n+]cc2)cc1", 2, True,
                 "two waves, spacing 0.42 V (benzyl-methyl / methyl viologen, tier A)"),
    "pyridinium": ("[n+;R;H0]", 1, False, "single reduction; parent reported irreversible"),
    "phenothiazine": ("c1cccc2c1[#16]c1ccccc1[#7]2", 1, False, "1e- oxidation (catholyte)"),
    "nitroxide": ("[#7;X3]([#6])([#6])[#8;X1;v1]", 1, False, "1e- oxidation (catholyte)"),
}
_PAT = {k: Chem.MolFromSmarts(v[0]) for k, v in FAMILIES.items() if v[0]}
_IMIDE = _PAT["mono-imide"]
_PYRIDINE_N = Chem.MolFromSmarts("[nX2;H0;+0]")    # Menshutkin-graftable pyridine N

# Merrifield attachment chemistry: nucleophile + ArCH2Cl -> ArCH2-Nu. Reactant atom [*:9]
# is the nucleophilic atom; the product gets the tether on it.
_T = "[CH2](c1ccc(C)cc1)"
GRAFTS = {
    "O-benzyl ether (phenol/enol, Williamson)":
        "[OX2H1:9]-[#6;!$([CX3]=O):2]>>" + _T + "[O:9]-[#6:2]",
    "N-benzyl (amine N-H)":
        "[NX3;!a;H1,H2;!$(N[#6]=O);!$(N[#16]=O);!$(N=*):9]>>" + _T + "[#7:9]",
    # aromatic or aliphatic imide N-H (RDKit perceives e.g. the PMDI imide N as aromatic [nH])
    "N-benzyl (imide N-H)":
        "[#7X3;H1:9]([#6:2]=[O:3])[#6:4]=[O:5]>>" + _T + "[#7:9]([#6:2]=[O:3])[#6:4]=[O:5]",
    "N-benzyl pyridinium (Menshutkin)":
        "[nX2;H0;+0:9]>>" + _T + "[n+:9]",
    "S-benzyl thioether":
        "[SX2H1:9]>>" + _T + "[S:9]",
    "benzyl ester (carboxylate)":
        "[CX3:2](=[O:3])[OX2H1:9]>>[CX3:2](=[O:3])[O:9]" + _T,
}
_RXN = {k: AllChem.ReactionFromSmarts(v) for k, v in GRAFTS.items()}
ALLOWED = {1, 5, 6, 7, 8, 9, 15, 16, 17, 35}   # H B C N O F P S Cl Br


def families(mol) -> list[str]:
    out = [k for k, p in _PAT.items() if mol.HasSubstructMatch(p)]
    if "mono-imide" in out and len(mol.GetSubstructMatches(_IMIDE)) >= 2:
        out = [f for f in out if f != "mono-imide"] + ["bis-imide"]
    if "viologen" in out:                       # a bipyridinium is not a lone pyridinium
        out = [f for f in out if f != "pyridinium"]
    return out


def graft_products(mol) -> list[tuple[str, str]]:
    """[(grafted canonical SMILES, chemistry)] — one per distinct attachment site."""
    seen, out = set(), []
    for chem, rxn in _RXN.items():
        for prods in rxn.RunReactants((mol,)):
            p = prods[0]
            try:
                Chem.SanitizeMol(p)
            except Exception:
                continue
            smi = Chem.MolToSmiles(p)
            if smi not in seen:
                seen.add(smi); out.append((smi, chem))
    return out


def _unassigned_stereo(mol) -> int:
    from rdkit.Chem import FindMolChiralCenters
    n = sum(1 for _, t in FindMolChiralCenters(mol, includeUnassigned=True,
                                                useLegacyImplementation=False) if t == "?")
    si = Chem.FindPotentialStereo(mol)
    n_db = sum(1 for s in si if s.type == Chem.StereoType.Bond_Double
               and s.specified == Chem.StereoSpecified.Unspecified)
    return n + n_db


def describe(grafted_smiles: str, n_e: int, faraday: float) -> dict:
    from redox.properties.capacity_and_proxies import _sa_prep, _sascorer
    m = Chem.MolFromSmiles(grafted_smiles)
    mw = Descriptors.MolWt(m)
    sa_mol, _, sa_note = _sa_prep(grafted_smiles)
    sa = _sascorer().calculateScore(sa_mol) if sa_mol is not None else None
    return dict(MW=round(mw, 3), n_e_nominal=n_e,
                capacity_nominal_mAh_g=round(n_e * faraday / (3.6 * mw), 1),
                SA_grafted=round(sa, 2) if sa is not None else None, sa_note=sa_note,
                charge_resting=Chem.GetFormalCharge(m), heavy_atoms=m.GetNumHeavyAtoms())


def identify(records, faraday: float, known_inchikeys=frozenset()) -> list[dict]:
    """records: iterable of dict(source_id, smiles, **extra). Returns one row per
    (molecule, graft site) for molecules in a FAMILIES family; rows carry the rule outcome
    (`status`) rather than being silently dropped."""
    rows = []
    for rec in records:
        m = Chem.MolFromSmiles(str(rec["smiles"]))
        if m is None:
            continue
        fams = families(m)
        if not fams and not m.HasSubstructMatch(_PYRIDINE_N):
            continue
        parent_key = Chem.MolToInchiKey(m)[:14]
        base = dict(rec, parent_smiles=Chem.MolToSmiles(m), parent_families="+".join(fams),
                    parent_in_library=parent_key in known_inchikeys)
        bad_el = sorted({a.GetSymbol() for a in m.GetAtoms() if a.GetAtomicNum() not in ALLOWED})
        grafts = graft_products(m)
        if not grafts:
            if fams:
                rows.append(dict(base, status="no graft handle", grafted_smiles=None,
                                 graft=None, families="+".join(fams),
                                 multi_e=any(FAMILIES[f][2] for f in fams)))
            continue
        for gsmi, chem in grafts:
            g = Chem.MolFromSmiles(gsmi)
            gf = families(g)                      # classify what is actually screened
            if not gf:
                continue
            multi = [f for f in gf if FAMILIES[f][2]]
            n_e = max(FAMILIES[f][1] for f in gf)
            # N-alkylating a ring N only builds a validated motif when it completes a
            # viologen (every registered Menshutkin graft does); elsewhere it turns the core
            # into an azinium cation (1e- pyridinium, or an unvalidated cationic quinone).
            menshutkin_off = chem.startswith("N-benzyl pyridinium") and "viologen" not in gf
            status = ("element not allowed: " + ",".join(bad_el) if bad_el else
                      "unassigned stereo" if _unassigned_stereo(g) else
                      "ring-N graft not forming a viologen" if menshutkin_off else
                      "1e- family" if not multi else "ok")
            rows.append(dict(base, grafted_smiles=gsmi, graft=chem, status=status,
                             families="+".join(gf), multi_e=bool(multi),
                             evidence="; ".join(FAMILIES[f][3] for f in gf),
                             grafted_in_library=Chem.MolToInchiKey(g)[:14] in known_inchikeys,
                             **describe(gsmi, n_e, faraday)))
    return rows


def precursor_from_frag(frag: str) -> str:
    """The un-grafted precursor of a config fragment: the [*:1] tether becomes H (or, on an
    azinium N, the neutral pyridine the Menshutkin alkylation starts from)."""
    m = Chem.RWMol(Chem.MolFromSmiles(frag))
    dummy = next(a for a in m.GetAtoms() if a.GetAtomicNum() == 0)
    nb = dummy.GetNeighbors()[0]
    m.RemoveAtom(dummy.GetIdx())
    if nb.GetIsAromatic() and nb.GetSymbol() == "N" and nb.GetFormalCharge() == 1:
        nb.SetFormalCharge(0)
    else:
        nb.SetNumExplicitHs(nb.GetNumExplicitHs() + 1)
    nb.SetNoImplicit(False)
    Chem.SanitizeMol(m)
    return Chem.MolToSmiles(m)


def rediscover(groups, faraday: float) -> list[dict]:
    """For each config group (id, frag): does applying the grafting rules to its precursor
    regenerate the grafted candidate? Returns rows with `rediscovered` True/False."""
    out = []
    for g in groups:
        target = Chem.MolToSmiles(Chem.molzip(Chem.MolFromSmiles(TETHER + "." + g["frag"])))
        pre = precursor_from_frag(g["frag"])
        prods = graft_products(Chem.MolFromSmiles(pre))
        hit = next((c for s, c in prods if s == target), None)
        fams = families(Chem.MolFromSmiles(target))       # classify the grafted species
        out.append(dict(id=g["id"], name=g.get("name"), precursor=pre, target=target,
                        rediscovered=hit is not None, graft=hit, n_sites=len(prods),
                        families="+".join(fams),
                        multi_e=any(FAMILIES[f][2] for f in fams)))
    return out
