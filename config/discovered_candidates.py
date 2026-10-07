"""Discovered candidates — the GROWING registry of grafted candidates found by the systematic
identification (src/redox/screening/discovery.py, scripts/mining/identify_candidates.py).

APPEND-ONLY: entries are added by scripts/mining/register_candidates.py and never rewritten,
so every computed result stays traceable to the entry it was computed for. Same schema as
config/starting_candidates.py (SCAFFOLD + frag with one [*:1] tether site; states =
(label, charge, mult, n_e), first = resting) plus `provenance` (source database + id, parent
SMILES, PubChem CID/title, graft chemistry, date) and an optional chemistry `flag`.

Names: the grafted species' real name derived from the PubChem parent name; as in every
candidate config, "benzyl" denotes the 4-methylbenzyl tether proxy of the polymer.

RUN protocol (identical for every entry; same as the starting-candidate batch):
  UMA (charge/spin) pre-opt -> DFT+SMD(MeCN) optimization r2SCAN-D4/def2-SVP(D) from
  RUN["nconf"] conformer seeds, each xtb-ALPB pre-optimized, lowest kept -> active-protocol
  single points wB97M-V/def2-TZVPD (redox.core.protocol.ACTIVE_SP) -> GFN2-xTB thermal ->
  Nelsen lambda cross points -> finalize chain.
"""

SCAFFOLD = "Cc1ccc(C[*:1])cc1"

RUN = dict(nconf=3, preopt="alpb")

GROUPS = [
    {'id': 'd3_05hfap',
     'name': '4-(benzylamino)-3-bromo-1,2-naphthoquinone',
     'family': 'quinone (n-type)',
     'frag': 'O=C1C(=O)c2ccccc2C(N[*:1])=C1Br',
     'states': [('neu', 0, 1, 0), ('red1', -1, 2, -1), ('red2', -2, 1, -2)],
     'provenance': {'source': 'D3TaLES',
                    'source_id': '05HFAP',
                    'parent_smiles': 'NC1=C(Br)C(=O)C(=O)c2ccccc21',
                    'pubchem_cid': 344278,
                    'pubchem_title': '4-Amino-3-bromonaphthalene-1,2-dione',
                    'graft': 'N-benzyl (amine N-H)',
                    'discovery_rule': 'redox.screening.discovery',
                    'registered': '2026-10-06'},
     'flag': 'stability: o-quinone (reactive Michael acceptor) with a ring C-Br substitution site; '
             'route as for aminoquinones'},
    {'id': 'd3_80izcq',
     'name': '4-{2-[4-(benzyloxy)phenyl]propan-2-yl}-1,2-benzoquinone',
     'family': 'quinone (n-type)',
     'frag': 'CC(C)(C1=CC(=O)C(=O)C=C1)c1ccc(O[*:1])cc1',
     'states': [('neu', 0, 1, 0), ('red1', -1, 2, -1), ('red2', -2, 1, -2)],
     'provenance': {'source': 'D3TaLES',
                    'source_id': '80IZCQ',
                    'parent_smiles': 'CC(C)(C1=CC(=O)C(=O)C=C1)c1ccc(O)cc1',
                    'pubchem_cid': 656690,
                    'pubchem_title': '4-(1-(4-Hydroxyphenyl)-1-methylethyl)-3,5-cyclohexadiene-1,2-dione',
                    'graft': 'O-benzyl ether (phenol/enol, Williamson)',
                    'discovery_rule': 'redox.screening.discovery',
                    'registered': '2026-10-06'},
     'flag': 'stability: unsubstituted o-benzoquinone positions (reactive Michael acceptor; '
             'o-quinones of this type are known to be unstable)'},
    {'id': 'd3_05stji',
     'name': '2-(benzylamino)-1,4-naphthoquinone',
     'family': 'quinone (n-type)',
     'frag': 'O=C1C=C(N[*:1])C(=O)c2ccccc21',
     'states': [('neu', 0, 1, 0), ('red1', -1, 2, -1), ('red2', -2, 1, -2)],
     'provenance': {'source': 'D3TaLES',
                    'source_id': '05STJI',
                    'parent_smiles': 'NC1=CC(=O)c2ccccc2C1=O',
                    'pubchem_cid': 72909,
                    'pubchem_title': '3-Aminonaphthoquinone',
                    'graft': 'N-benzyl (amine N-H)',
                    'discovery_rule': 'redox.screening.discovery',
                    'registered': '2026-10-06'},
     'flag': 'route: the 2-amino N is conjugated to the quinone (poor nucleophile); the same '
             'product is normally made by adding a benzylamine / aminomethyl resin to the quinone, '
             'not by N-benzylation'},
    {'id': 'd3_05clzk',
     'name': '2-(benzyloxy)-3-chloro-1,4-naphthoquinone',
     'family': 'quinone (n-type)',
     'frag': 'O=C1C(Cl)=C(O[*:1])C(=O)c2ccccc21',
     'states': [('neu', 0, 1, 0), ('red1', -1, 2, -1), ('red2', -2, 1, -2)],
     'provenance': {'source': 'D3TaLES',
                    'source_id': '05CLZK',
                    'parent_smiles': 'O=C1C(O)=C(Cl)C(=O)c2ccccc21',
                    'pubchem_cid': 73711,
                    'pubchem_title': '2-Chloro-3-hydroxy-1,4-naphthoquinone',
                    'graft': 'O-benzyl ether (phenol/enol, Williamson)',
                    'discovery_rule': 'redox.screening.discovery',
                    'registered': '2026-10-06'},
     'flag': 'stability: the ring C-Cl is a nucleophilic vinylic substitution site (amines, '
             'thiols, alkoxides) - may not survive grafting or cycling'},
    {'id': 'd3_05hdcy',
     'name': '8-(benzyloxy)-2-chloro-7-ethyl-6-methoxy-1,4-naphthoquinone',
     'family': 'quinone (n-type)',
     'frag': 'CCc1c(OC)cc2c(c1O[*:1])C(=O)C(Cl)=CC2=O',
     'states': [('neu', 0, 1, 0), ('red1', -1, 2, -1), ('red2', -2, 1, -2)],
     'provenance': {'source': 'D3TaLES',
                    'source_id': '05HDCY',
                    'parent_smiles': 'CCc1c(OC)cc2c(c1O)C(=O)C(Cl)=CC2=O',
                    'pubchem_cid': 71406861,
                    'pubchem_title': '2-Chloro-7-ethyl-8-hydroxy-6-methoxynaphthalene-1,4-dione',
                    'graft': 'O-benzyl ether (phenol/enol, Williamson)',
                    'discovery_rule': 'redox.screening.discovery',
                    'registered': '2026-10-06'},
     'flag': 'stability: ring C-Cl substitution site; grafting the peri 8-OH removes its '
             'intramolecular H-bond to C1=O (shifts potentials)'},
    {'id': 'd3_80mnnp',
     'name': '2-(benzyloxy)-3-bromo-1,4-naphthoquinone',
     'family': 'quinone (n-type)',
     'frag': 'O=C1C(Br)=C(O[*:1])C(=O)c2ccccc21',
     'states': [('neu', 0, 1, 0), ('red1', -1, 2, -1), ('red2', -2, 1, -2)],
     'provenance': {'source': 'D3TaLES',
                    'source_id': '80MNNP',
                    'parent_smiles': 'O=C1C(O)=C(Br)C(=O)c2ccccc21',
                    'pubchem_cid': 344276,
                    'pubchem_title': '2-Bromo-3-hydroxynaphthalene-1,4-dione',
                    'graft': 'O-benzyl ether (phenol/enol, Williamson)',
                    'discovery_rule': 'redox.screening.discovery',
                    'registered': '2026-10-06'},
     'flag': 'stability: the ring C-Br is a nucleophilic vinylic substitution site (amines, '
             'thiols, alkoxides) - may not survive grafting or cycling'},
    {'id': 'd3_80txmo',
     'name': '2-[benzyl(ethyl)amino]-1,4-naphthoquinone',
     'family': 'quinone (n-type)',
     'frag': 'CCN(C1=CC(=O)c2ccccc2C1=O)[*:1]',
     'states': [('neu', 0, 1, 0), ('red1', -1, 2, -1), ('red2', -2, 1, -2)],
     'provenance': {'source': 'D3TaLES',
                    'source_id': '80TXMO',
                    'parent_smiles': 'CCNC1=CC(=O)c2ccccc2C1=O',
                    'pubchem_cid': 312577,
                    'pubchem_title': '2-(Ethylamino)naphthalene-1,4-dione',
                    'graft': 'N-benzyl (amine N-H)',
                    'discovery_rule': 'redox.screening.discovery',
                    'registered': '2026-10-06'},
     'flag': 'route: tertiary enaminone; sterically hindered N-benzylation of a conjugated (weakly '
             'nucleophilic) amine'},
    {'id': 'd3_80fsak',
     'name': '2-(benzyloxy)-3-methyl-1,4-naphthoquinone (O-benzyl phthiocol)',
     'family': 'quinone (n-type)',
     'frag': 'CC1=C(O[*:1])C(=O)c2ccccc2C1=O',
     'states': [('neu', 0, 1, 0), ('red1', -1, 2, -1), ('red2', -2, 1, -2)],
     'provenance': {'source': 'D3TaLES',
                    'source_id': '80FSAK',
                    'parent_smiles': 'CC1=C(O)C(=O)c2ccccc2C1=O',
                    'pubchem_cid': 10221,
                    'pubchem_title': 'Phthiocol',
                    'graft': 'O-benzyl ether (phenol/enol, Williamson)',
                    'discovery_rule': 'redox.screening.discovery',
                    'registered': '2026-10-06'},
     'flag': 'route: O- vs C-alkylation selectivity of the acidic 2-hydroxy-1,4-naphthoquinone '
             'enol (as for lawsone)'},
# >>> append new entries above this line (register_candidates.py) <<<
]
