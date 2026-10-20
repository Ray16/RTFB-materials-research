#!/usr/bin/env python
"""Independent structural admission for stationary polymer electronic records."""
import copy
import hashlib
import json
import math
from pathlib import Path

from redox.core.common import ROOT,EV_KJ
from redox.core.io import atomic_write_json
from redox.core.protocol import load_record
from redox.build.geometry import geometry_stereo_status
from redox.build.polymer_contacts import contact_metrics
from redox.validation.protonation import hydrogen_topology_status
from ase.io import read


def main():
    out=ROOT/'results/polymer_models'
    source=out/'interaction_results.json';feed=json.loads(source.read_text());checked=copy.deepcopy(feed)
    admitted={};held=[]
    for row in checked['models']:
        if row['final_xyz'] is None:continue
        path=Path(row['final_xyz']);sd=path.parent;issues=[]
        stereo=geometry_stereo_status(row['smiles'],path)
        if stereo not in ('ok','not_applicable'):issues.append('stereochemistry:'+stereo)
        topology=hydrogen_topology_status(row['smiles'],path)
        if topology!='ok':issues.append('connectivity:'+topology)
        raw=json.loads((sd/'result.json').read_text())
        if raw['charge']!=0 or raw['mult']!=1:issues.append('charge_or_spin')
        record=load_record(sd,0,1)
        if record is None or not math.isclose(record['e_smd_eV'],row['E_smd_eV'],rel_tol=0,abs_tol=1e-9):issues.append('active_SP')
        contacts=contact_metrics(read(path).positions,row['pi_groups'])
        expected=2 if row['model_type']=='two_nonbonded_segments' else row['repeats']
        if len(row['pi_groups'])!=expected:issues.append('redox_core_count')
        row['structural_review']=dict(stereo=stereo,topology=topology,issues=issues)
        if issues:
            held.append(dict(key=row['key'],issues=issues));row['reported_E_smd_eV']=row['E_smd_eV'];row['E_smd_eV']=None;row['stage']='held_structural_review'
        else:
            row['contact_metrics']=contacts;row['stacked_pairs']=sum(p['stack_geometry'] for p in contacts)
            admitted[row['key']]=row
    by_model={}
    for row in checked['models']:by_model.setdefault(row['model_key'],{})[row['target']]=row
    for comparison in checked['folding']:
        states=by_model[comparison['model_key']];a,b=states['open'],states['stacked']
        valid=a['key'] in admitted and b['key'] in admitted and a['stacked_pairs']==0 and b['stacked_pairs']>0
        comparison['stacked_minus_open_deltaE_kJ_mol']=(b['E_smd_eV']-a['E_smd_eV'])*EV_KJ if valid else None
        comparison['structurally_reviewed']=valid
    monomers={r['candidate_id']:r for r in checked['models'] if r['target']=='monomer'}
    for comparison in checked['association']:
        row=next(r for r in checked['models'] if r['model_key']==comparison['model_key'])
        mono=monomers[row['candidate_id']]
        valid=row['key'] in admitted and mono['key'] in admitted
        comparison['association_deltaE_kJ_mol']=(row['E_smd_eV']-2*mono['E_smd_eV'])*EV_KJ if valid else None
        comparison['structurally_reviewed']=valid
    checked['status'].update(structurally_reviewed_DFT_cases=len(admitted),stationary_electronic_cases=len(admitted),
        structural_holds=held,sampled_stack_open_comparisons=sum(r['stacked_minus_open_deltaE_kJ_mol'] is not None for r in checked['folding']),
        electronic_association_comparisons=sum(r['association_deltaE_kJ_mol'] is not None for r in checked['association']))
    campaign=ROOT/checked['status']['campaign']
    plan=json.loads((campaign/'plan.json').read_text())
    pilot_ready=(campaign/'live_inputs'/plan['pilot_key']/'task.json').exists()
    if not pilot_ready and not (campaign/'pilot/status.json').exists():
        checked['status']['memory_pilot']='awaiting_largest_trimer_warm_start'
    checked['status']['prepared_DFT_inputs']=len(list((campaign/'live_inputs').glob('*/task.json')))
    checked['review']=dict(raw_source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        review_source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        scope='Connectivity, backbone stereochemistry, declared charge/spin, active SP and independently recomputed geometric contacts; not Hessian minima or free energies')
    atomic_write_json(out/'interaction_checked.json',checked)
    atomic_write_json(out/'interaction_status.json',checked['status'])
    atomic_write_json(out/'interaction_review.json',dict(**checked['review'],admitted=len(admitted),held=held))
    print('Structural review:',len(admitted),'stationary electronic records;',len(held),'held')


if __name__=='__main__':main()
