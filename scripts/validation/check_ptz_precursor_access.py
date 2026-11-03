"""Check exact precursor-listing identity without promoting an attachment route."""
import argparse
import json
import os
from pathlib import Path

from rdkit import Chem

from redox.core.common import ROOT
from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json
from redox.workflows.molecular_native_namespace import load


def matching_precursors(panel, properties, supplier):
    compound=properties['PropertyTable']['Properties'][0]
    sid=supplier['InformationList']['Information'][0]
    if sid['SID']!=468592587 or sid['CID']!=[compound['CID']] or compound['CID']!=98841:
        raise ValueError('The primary supplier deposited substance must identify this exact compound')
    def canonical(smiles):
        mol=Chem.MolFromSmiles(smiles)
        if mol is None:raise ValueError('Valid exact precursor graph required')
        return Chem.MolToSmiles(mol,isomericSmiles=True)
    precursor=canonical(compound['SMILES'])
    if Chem.MolToInchiKey(Chem.MolFromSmiles(precursor))!=compound['InChIKey']:
        raise ValueError('Deposited precursor graph and structure identifier differ')
    return [row for row in panel if canonical(row['precursor_smiles'])==precursor]


def audit(folder,previous,current):
    def read(name):return json.loads((folder/name).read_text())
    source=read('source_review.json')
    if any(sha256_file(p)!=h for p,h in source['source_sha256'].items()):
        raise ValueError('Primary precursor or analogue evidence changed')
    for name,original in [('precursor_properties.json','PubChem_98841_properties.json'),
            ('supplier_sid_cids.json','TCI_PubChem_substance_CIDs.json'),('primary_doi_metadata.json','Cao2024_crossref.json')]:
        bound=[h for p,h in source['source_sha256'].items() if Path(p).name==original]
        if len(bound)!=1 or sha256_file(folder/name)!=bound[0]:
            raise ValueError('Curated precursor/metadata copy differs from the consulted primary source')
    old=json.loads((previous/'registry.json').read_text())['candidates']
    new=json.loads((current/'registry.json').read_text())['candidates']
    if sha256_file(folder/'route_evidence_before.json')!=sha256_file(previous/'inputs/route_evidence.json'):
        raise ValueError('The earlier route catalog must match the preserved previous pool')
    oldmap={r['graph_key']:r for r in old};newmap={r['graph_key']:r for r in new}
    if len(old)!=912 or len(new)!=912 or oldmap.keys()!=newmap.keys():
        raise ValueError('All 912 existing exact graphs must be retained')
    for key,row in oldmap.items():
        if {k:v for k,v in row.items() if k!='route_decisions'}!={k:v for k,v in newmap[key].items() if k!='route_decisions'}:
            raise ValueError('Precursor evidence cannot change graphs, aliases, evaluations or material decisions')
    panel=json.loads((current/'inputs/ptz_direct_ring_2026-11-02.json').read_text())
    bound=[h for p,h in source['source_sha256'].items() if Path(p).name=='ptz_direct_ring_2026-11-02.json']
    if len(bound)!=1 or sha256_file(current/'inputs/ptz_direct_ring_2026-11-02.json')!=bound[0]:
        raise ValueError('Exact design panel changed after precursor matching')
    matches=matching_precursors(panel,read('precursor_properties.json'),read('supplier_sid_cids.json'))
    if len(matches)!=1 or matches[0]['id']!=source['candidate_alias']:
        raise ValueError('Only the exact 2-cyano N-H precursor may receive this listing evidence')
    row=matches[0];evidence=json.loads((current/'inputs/route_evidence.json').read_text())['candidates']
    gate=load(ROOT/'scripts/pipeline/screen_grafting_routes.py','independent_ptz_access_gate')
    decision=gate.screen(row,evidence,{})
    if (gate.exact_route_binding(row)!=source['route_binding'] or evidence[row['id']]['route_binding']!=source['route_binding']
            or evidence[row['id']]['exact_precursor_access_reviewed'] is not True
            or decision['route_status']!='hold' or decision['proceed_to_full_redox']
            or 'exact_precursor_access_reviewed' in decision['missing_route_evidence']
            or set(decision['missing_route_evidence'])!={'attachment_selectivity_reviewed','balanced_route_reviewed','redox_core_retention_reviewed'}):
        raise ValueError('Exact bond binding and the remaining route hold must be retained')
    if {k:v for k,v in evidence.items() if k!=row['id']}!={k:v for k,v in read('route_evidence_before.json')['candidates'].items() if k!=row['id']}:
        raise ValueError('Unrelated route evidence changed')
    summary=json.loads((current/'summary.json').read_text())
    if (summary['route_counts']!={'proceed':3,'hold':909} or summary['full_material_admissions']!=0
            or summary['new_quantum_jobs']!=0 or any(r['material_decision']!='pending' for r in new)):
        raise ValueError('A precursor listing supplies no route/material or property admission')
    return dict(status='exact_precursor_access_and_registry_preservation_checked',candidate_alias=row['id'],
        exact_graphs_preserved=912,prior_evaluations_unchanged=True,new_quantum_calculations=0,
        route_counts=summary['route_counts'],full_material_admissions=0,remaining_route_gaps=decision['missing_route_evidence'],
        precursor_access_reviewed=True,full_redox_campaign_admitted=False,properties_transferred=False,
        source_review_sha256=sha256_file(folder/'source_review.json'),registry_sha256=sha256_file(current/'registry.json'),
        producer_sha256=sha256_file(__file__))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('review','previous','current','out'):p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args()
    from redox.workflows.resources import allocation_guard,host_capacity,cpu_allocation
    with allocation_guard():host_capacity(1,minimum_ram_gib=4);os.sched_setaffinity(0,cpu_allocation(1))
    if a.out.exists():raise FileExistsError('Preserve earlier precursor access reviews')
    result=audit(a.review.resolve(),a.previous.resolve(),a.current.resolve());atomic_write_json(a.out,result)
    print(json.dumps(result))
