"""Replay the exact registry extension without promoting graph/route QC to properties."""
from collections import Counter
import hashlib
import importlib.util
import json
from pathlib import Path

from rdkit import Chem

ROOT=Path(__file__).resolve().parents[1]


def reconstruct(folder):
    def read(name):return json.loads((folder/('registry_update_'+name+'.json')).read_text())
    old=read('old')['candidates'];new=read('current')['candidates'];summary=read('summary')
    preservation=read('preservation');panel=read('panel');evidence=read('route_evidence')['candidates']
    previous={row['graph_key']:row for row in old};current={row['graph_key']:row for row in new}
    if len(previous)!=len(old) or len(current)!=len(new) or not previous.keys()<=current.keys():
        raise ValueError('Exact graph identities must remain unique and all prior graphs preserved')
    for key,row in current.items():
        mol=Chem.MolFromSmiles(row['product_smiles'])
        if (mol is None or Chem.MolToSmiles(mol,isomericSmiles=True)!=row['product_smiles']
                or hashlib.sha256(row['product_smiles'].encode()).hexdigest()!=key
                or row['material_decision']!='pending'):
            raise ValueError('Exact graph or field-specific material scope differs')
    for key,row in previous.items():
        if (row['product_smiles']!=current[key]['product_smiles'] or row['evaluations']!=current[key]['evaluations']
                or not set(row['aliases'])<=set(current[key]['aliases'])):
            raise ValueError('Registry expansion changed a prior graph, alias or reported evaluation')
    added=[r for key,r in current.items() if key not in previous]
    if (len(old)!=904 or len(new)!=912 or len(added)!=8 or len(panel)!=9
            or preservation['old_graphs_preserved']!=len(old) or preservation['new_graphs']!=len(added)
            or summary['exact_graphs']!=len(new) or summary['route_counts']!=dict(Counter(r['route_status'] for r in new))
            or summary['route_counts']!={'hold':909,'proceed':3}
            or summary['full_material_admissions']!=0 or summary['new_quantum_jobs']!=0
            or summary['raw_quantum_records_revalidated']):
        raise ValueError('Exact extension census and graph/route-only scope required')
    if any(r['evaluations'] or r['route_status']!='hold' for r in added):
        raise ValueError('No measured/computed property or route pass transfers to the new designs')
    module=importlib.util.spec_from_file_location('frozen_registry_route_gate',ROOT/'scripts/pipeline/screen_grafting_routes.py')
    gate=importlib.util.module_from_spec(module);module.loader.exec_module(gate)
    if sum(len(r['recipes']) for r in panel)!=17:raise ValueError('Retain all symmetry-equivalent generation recipes')
    for row in panel:
        if (row['properties_computed'] or row['route_binding']!=gate.exact_route_binding(row)
                or row['states']!=[['neu',0,1,0],['ox',1,2,1]]):
            raise ValueError('Exact N-linked pairs and one-electron state hypotheses required')
        decision=gate.screen(row,evidence,{'candidates':{}})
        if decision['proceed_to_full_redox']!=(row['id']=='phenothiazine'):
            raise ValueError('Only the reviewed unsubstituted baseline may pass this panel')
    baseline=next(r for r in current.values() if 'phenothiazine' in r['aliases'])
    if baseline['route_status']!='proceed' or baseline['evaluations']:
        raise ValueError('The route-qualified baseline still has no admitted exact-graft properties here')
    return dict(exact_graphs=len(new),prior_graphs_preserved=len(old),new_uncomputed_graphs=len(added),
        route_passes=3,route_holds=909,full_material_admissions=0,prior_evaluations_unchanged=True,
        raw_quantum_records_revalidated=False,new_quantum_calculations=False)


def values(data):
    return ''.join('\\newcommand{\\'+name+'}{'+str(value)+'}\n' for name,value in [
        ('UnifiedRegistryCount',data['exact_graphs']),('NewDirectPTZDesigns',data['new_uncomputed_graphs']),
        ('ReviewedRouteCount',data['route_passes'])])


def audit(folder):
    result=reconstruct(folder)
    if (folder/'registry_update_values.tex').read_text()!=values(result):
        raise ValueError('Registry count macros differ from their frozen graphs and route decisions')
    return result
