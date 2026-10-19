"""Independently verify the identities and electron counting of all state drawings."""
import hashlib
import json
from pathlib import Path

import networkx as nx
from rdkit import Chem
from pypdf import PdfReader

ROOT=Path(__file__).resolve().parents[4]
OUT=Path(__file__).resolve().parent
data=json.loads((ROOT/'results/candidates/figures/data/redox_states.json').read_text())
proposal=json.loads((ROOT/'results/discovery/vinyl_monomers/proposals.json').read_text())
proposals={r['id']:r for r in proposal['records']}
checks=[]


def check(name,passed):checks.append(dict(check=name,passed=bool(passed)))


def graph(mol):
    g=nx.Graph()
    for atom in mol.GetAtoms():g.add_node(atom.GetIdx(),element=atom.GetAtomicNum(),hydrogens=atom.GetTotalNumHs(),isotope=atom.GetIsotope())
    for bond in mol.GetBonds():g.add_edge(bond.GetBeginAtomIdx(),bond.GetEndAtomIdx())
    return g


for state in data['records']:
    prefix=state['id']+'/'+state['state']+' '
    source=Chem.MolFromSmiles(state['source_smiles']);drawing=Chem.MolFromSmiles(state['depicted_smiles'])
    same=nx.is_isomorphic(graph(source),graph(drawing),node_match=nx.algorithms.isomorphism.categorical_node_match(['element','hydrogens','isotope'],[None,None,None]))
    check(prefix+'independent atom/H/connectivity isomorphism',same)
    charge=sum(atom.GetFormalCharge() for atom in drawing.GetAtoms())
    radicals=sum(atom.GetNumRadicalElectrons() for atom in drawing.GetAtoms())
    nuclei=sum(atom.GetAtomicNum() for atom in Chem.AddHs(drawing).GetAtoms())
    source_electrons=sum(atom.GetAtomicNum() for atom in Chem.AddHs(source).GetAtoms())-Chem.GetFormalCharge(source)
    check(prefix+'molecular charge',charge==state['charge'])
    check(prefix+'electron gain without proton transfer',nuclei-charge-source_electrons==-state['n_e'])
    check(prefix+'radical count and multiplicity parity',radicals==state['radicals']==state['mult']-1 and (nuclei-charge-radicals)%2==0)
    if state['scope']=='computed_selected_spin_ladder':
        paths=state['evidence_sha256'];check(prefix+'computed source evidence exists',len(paths)==1)
        for p,digest in paths.items():
            raw=json.loads(Path(p).read_text());check(prefix+'matching computed state',raw['charge']==charge and raw['mult']==state['mult'])
            check(prefix+'unchanged source bytes',hashlib.sha256(Path(p).read_bytes()).hexdigest()==digest)
    else:
        declared=next(s for s in proposals[state['id']]['states'] if s[0]==state['state'])
        check(prefix+'proposed label, no computed evidence',declared[1]==charge and declared[2]==state['mult'] and not state['evidence_sha256'])
for name in ('redox_states_computed.pdf','redox_states_vinyl.pdf'):
    page=PdfReader(ROOT/'results/candidates/figures/structures'/name).pages[0]
    check(name+' vector export',not list(page.images))
    check(name+' publication width',abs(float(page.mediabox.width)-7.2*72)<.05)
check('Combined two-page state comparison',len(PdfReader(ROOT/'results/candidates/figures/structures/redox_states.pdf').pages)==2)
check('Nine computed and nine proposed state depictions',sum(r['scope']=='computed_selected_spin_ladder' for r in data['records'])==9 and sum(r['scope']=='proposed_uncomputed_ladder' for r in data['records'])==9)
result=dict(status='independently_checked_formal_state_schematics',states=len(data['records']),checks=len(checks),failed_checks=sum(not c['passed'] for c in checks),
            calculated_spin_density=False,vinyl_quantum_properties_computed=False,details=checks)
(OUT/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k!='details'}))
assert not result['failed_checks']
