#!/usr/bin/env python
"""Reconstruct a complete matched neighboring-graft cycle from native states.

Shared base/precursor terms cancel from nonadditivity. These finite-minimum
cycles do not isolate motif binding, electronic coupling or an occupancy law.
"""
import argparse
import json
import math
import os
from pathlib import Path
import time

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api

EV_KJ = 96.48533212331002


def cycle(values, constant):
    if set(values) != {'empty','site2','site3','pair'} or not all(math.isfinite(v) for v in [*values.values(),constant]):
        raise ValueError('Four finite matched state free energies and a finite reaction constant are required')
    g0,g2,g3,g23 = (values[k] for k in ('empty','site2','site3','pair'))
    nonadditivity = math.fsum([g23,g0,-g2,-g3])*EV_KJ
    result = dict(graft_2_on_empty_kJ_mol=math.fsum([g2,-g0,constant])*EV_KJ,
        graft_3_on_empty_kJ_mol=math.fsum([g3,-g0,constant])*EV_KJ,
        graft_2_given_3_kJ_mol=math.fsum([g23,-g3,constant])*EV_KJ,
        graft_3_given_2_kJ_mol=math.fsum([g23,-g2,constant])*EV_KJ,
        direct_double_graft_kJ_mol=math.fsum([g23,-g0,2*constant])*EV_KJ,
        deltaDeltaG_kJ_mol=nonadditivity)
    differences=[result['graft_2_given_3_kJ_mol']-result['graft_2_on_empty_kJ_mol'],
                 result['graft_3_given_2_kJ_mol']-result['graft_3_on_empty_kJ_mol']]
    paths=[result['graft_2_on_empty_kJ_mol']+result['graft_3_given_2_kJ_mol'],
           result['graft_3_on_empty_kJ_mol']+result['graft_2_given_3_kJ_mol']]
    if (max(abs(v-nonadditivity) for v in differences)>1e-7
            or max(abs(v-result['direct_double_graft_kJ_mol']) for v in paths)>1e-7):
        raise ValueError('Conditional and direct thermodynamic paths do not close')
    return result


def analyze(first_path, pair_campaign, pair_review, out):
    first=json.loads(first_path.read_text());review=json.loads(pair_review.read_text())
    api.verify_sources(first['source_sha256']);api.verify_sources(review['source_sha256'])
    if (first['status']!='available_first_site_grafting_independently_reconstructed'
            or review['status']!='available_exact_polymer_native_components_independently_reconstructed'
            or review['checked_components']!=1 or len(review['states'])!=1
            or not review['normal_admission_rechecked'] or review['QM_recomputed']):
        raise ValueError('Checked first-site and completed paired-native component evidence required')
    first_api=api.load_module('matched_neighbor_independent_state_reader',ROOT/'scripts/validation/analyze_first_site_grafting.py')
    admitted,invalid=api.admitted_records(pair_campaign,'thermal')
    if invalid or len(admitted)!=1:raise ValueError('Exactly one normally admitted paired native state is required')
    record=next(iter(admitted.values()));task=record['task'];raw=review['states'][0]
    if task['occupied_sites']!=[2,3] or Path(record['state_dir'])!=Path(raw['state_dir']):
        raise ValueError('Paired state must be the exact reviewed sites-2/3 model')
    pair=first_api.state_values(Path(record['state_dir']),0,1,raw['G_eV'])
    if any(abs(pair['G_by_cutoff_eV'][str(k)]-raw['G_by_rotor_cutoff_eV'][str(k)])>1e-8 for k in (25,50,100)):
        raise ValueError('Paired cutoff partitions differ from their independent reconstruction')
    sources={**first['source_sha256'],**review['source_sha256'],**pair['source_sha256']}
    states={'pair':dict(pair,smiles=task['smiles'],task=task)}
    empty=first['empty_reference'];states['empty']=dict(first_api.state_values(Path(empty['state_dir']),0,1,empty['G_eV']),smiles=empty['smiles'])
    for single in first['singles']:
        if (single['task']['candidate_id']!=task['candidate_id']
                or single['task']['backbone_stereo_fingerprint']!=task['backbone_stereo_fingerprint']
                or single['task']['reactive_sites']!=task['reactive_sites']):
            raise ValueError('Candidate/backbone configuration/reactive-site identity differs')
        name='site'+str(single['site'])
        states[name]=dict(first_api.state_values(Path(single['state_dir']),0,1,single['G_eV']),
            smiles=single['task']['smiles'],task=single['task'])
    if set(states)!={'empty','site2','site3','pair'}:raise ValueError('Four matched occupancy graphs are required')
    balances={'nonadditivity':first_api.balance([states['pair']['smiles'],states['empty']['smiles']],
        [states['site2']['smiles'],states['site3']['smiles']])}
    reagents=first['reagent_states']
    for site,other in [(2,3),(3,2)]:
        left=[states['site'+str(other)]['smiles'],reagents['precursor']['smiles'],reagents['base']['smiles']]
        balances['graft'+str(site)+'_given'+str(other)+'_free_ions']=first_api.balance(left,
            [states['pair']['smiles'],reagents['protonated']['smiles'],reagents['chloride']['smiles']])
        balances['graft'+str(site)+'_given'+str(other)+'_contact_pair']=first_api.balance(left,
            [states['pair']['smiles'],reagents['pair']['smiles']])
    comparisons={}
    for cutoff in (25,50,100):
        key=str(cutoff);values={name:r['G_by_cutoff_eV'][key] for name,r in states.items()}
        models={name:cycle(values,constant) for name,constant in first['shared_constant_eV_by_cutoff'][key].items()}
        if abs(models['free_ions']['deltaDeltaG_kJ_mol']-models['contact_pair']['deltaDeltaG_kJ_mol'])>1e-10:
            raise ValueError('The common base/byproduct model must cancel from nonadditivity')
        for single in first['singles']:
            for name in models:
                old=single['by_cutoff'][key]['first_graft_kJ_mol'][name]
                if abs(models[name]['graft_'+str(single['site'])+'_on_empty_kJ_mol']-old)>1e-7:
                    raise ValueError('New cycle does not retain the independently checked first-site reaction')
        comparisons[key]=dict(state_G_eV=values,models=models,
            deltaDeltaG_kJ_mol=models['free_ions']['deltaDeltaG_kJ_mol'])
    electronic=math.fsum([states['pair']['E_eV'],states['empty']['E_eV'],-states['site2']['E_eV'],-states['site3']['E_eV']])*EV_KJ
    if abs(electronic-first['electronic_neighbor_cycle_kJ_mol'])>1e-7:
        raise ValueError('The electronic reference cycle changed')
    for row in states.values():sources.update(row['source_sha256'])
    for p in (Path(__file__).resolve(),first_path,pair_review,pair_campaign/'plan.json'):
        sources[str(p)]=api.digest(p)
    api.verify_sources(sources)
    for row in states.values():row.pop('source_sha256')
    report=dict(status='exact_candidate_neighbor_grafting_cycle_independently_reconstructed',
        checked=time.time(),candidate_id=first['candidate_id'],solvent='acetonitrile',temperature_K=298.15,
        base='DBU',conversion_kJ_mol_per_eV=EV_KJ,primary_cutoff_cm=50,
        states=states,by_cutoff=comparisons,electronic_deltaDeltaE_kJ_mol=electronic,
        deltaDeltaG_kJ_mol=comparisons['50']['deltaDeltaG_kJ_mol'],balances=balances,
        shared_reagent_constants_eV=first['shared_constant_eV_by_cutoff'],
        paired_component_QC=raw,empty_reference_selection=empty['selection'],
        conformers_per_graph=1,conformer_convergence_established=False,synthesis_admission=False,
        isolated_motif_interaction_computed=False,electronic_coupling_computed=False,
        graft_loading_distribution_inferred=False,source_sha256=sources,
        scope='Finite matched neutral-chain conditional grafting nonadditivity; whole backbone/solvent response included. Base/precursor terms cancel, absolute DBU route retains speciation qualifications.')
    out.mkdir(parents=True,exist_ok=False);atomic_write_json(out/'summary.json',report)
    print(json.dumps(dict(deltaDeltaG_by_cutoff_kJ_mol={k:v['deltaDeltaG_kJ_mol'] for k,v in comparisons.items()},
        primary_DBU_cycles=comparisons['50']['models'],conformer_convergence_established=False)))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ('first-site','pair-campaign','pair-review','out'):parser.add_argument('--'+name,type=Path,required=True)
    args=parser.parse_args()
    from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity
    with allocation_guard():host_capacity(1);os.sched_setaffinity(0,cpu_allocation(1))
    analyze(args.first_site.resolve(),args.pair_campaign.resolve(),args.pair_review.resolve(),args.out.resolve())


if __name__=='__main__':main()
