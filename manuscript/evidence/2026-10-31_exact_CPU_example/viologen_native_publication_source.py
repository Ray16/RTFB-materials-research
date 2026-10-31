#!/usr/bin/env python3
"""Assemble exact molecular redox states with checked Fc and explicit missing axes."""
import argparse
import importlib.util
import json
import math
from pathlib import Path
import shutil

from redox.core.common import ROOT,EV_KJ
from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json
from redox.core.protocol import ACTIVE_SP
from redox.qm.native_thermal import NATIVE_THERMAL
from redox.screening.feasibility import assess_candidate
from redox.validation.archived_sources import ArchivedSources
from redox.validation.molecular_native import MolecularSources


def ladder(states,ferrocene):
    """Reference exact one-electron charge steps; no fitted offset or imputation."""
    if len(states)<2 or len({r['charge'] for r in states})!=len(states):
        raise ValueError('One admitted state per distinct charge required')
    if any(type(r['charge']) is not int or type(r['mult']) is not int or r['mult']<1 for r in states):
        raise ValueError('Explicit integral charge and positive multiplicity required')
    if len({r['product_inchikey'] for r in states})!=1:
        raise ValueError('All redox states must have the same exact product identity')
    cutoffs={'25','50','100'}
    rows=[*states,*ferrocene.values()]
    if (set(ferrocene)!={'neu','ox'} or any(set(r['G_by_cutoff_eV'])!=cutoffs for r in rows)
            or any(type(v) not in (float,int) or not math.isfinite(v) for r in rows for v in r['G_by_cutoff_eV'].values())):
        raise ValueError('Finite native free energies at all three cutoffs required')
    ordered=sorted(states,key=lambda r:r['charge'],reverse=True);events=[]
    for oxidized,reduced in zip(ordered,ordered[1:]):
        if oxidized['charge']-reduced['charge']!=1:raise ValueError('One-electron adjacent charge steps required')
        values={c:math.fsum([oxidized['G_by_cutoff_eV'][c],-reduced['G_by_cutoff_eV'][c],
                 -ferrocene['ox']['G_by_cutoff_eV'][c],ferrocene['neu']['G_by_cutoff_eV'][c]]) for c in sorted(cutoffs,key=int)}
        events.append(dict(event=oxidized['state']+'->'+reduced['state'],electrons=1,
            initial_charge=oxidized['charge'],final_charge=reduced['charge'],
            initial_multiplicity=oxidized['mult'],final_multiplicity=reduced['mult'],
            potential_V_vs_Fc_by_cutoff=values))
    result=dict(events=events,fitted_offset_V=0.,nominal_temperature_K=298.15,solvent='acetonitrile',
        declared_spin_models=True,experimental_exact_graft_accuracy_validated=False)
    if len(states)==3:
        spacing={c:events[0]['potential_V_vs_Fc_by_cutoff'][c]-events[1]['potential_V_vs_Fc_by_cutoff'][c] for c in cutoffs}
        disprop={c:math.fsum([ordered[0]['G_by_cutoff_eV'][c],ordered[2]['G_by_cutoff_eV'][c],-2*ordered[1]['G_by_cutoff_eV'][c]]) for c in cutoffs}
        if any(abs(spacing[c]-disprop[c])>1e-8 for c in cutoffs):raise ValueError('Spacing and balanced disproportionation differ')
        result.update(spacing_V_by_cutoff=spacing,disproportionation_kJ_mol_by_cutoff={c:v*EV_KJ for c,v in disprop.items()})
    return result


def publish(status_path,reference_path,archive,registry_path,policy_path,output,electrode_pool):
    if output.exists():raise FileExistsError('Preserve earlier exact-candidate publications')
    status=json.loads(status_path.read_text());reviews=[];numeric=[]
    if (len(status['entries'])!=status['expected'] or status['checked']!=status['expected']
            or status['terminal']!=status['expected']):
        raise ValueError('The complete requested native-state set must pass independent review')
    for row in status['entries']:
        if row['status']!='exact_molecular_native_G_independently_checked':raise ValueError('Every requested state must pass independent review')
        path=Path(row['review'])
        if sha256_file(path)!=row['review_sha256']:raise ValueError('State review changed')
        review=json.loads(path.read_text());spec=json.loads(Path(review['namespace']).read_text());reader=MolecularSources(spec)
        if (review['status']!='exact_molecular_native_G_independently_checked'
                or review['resource_pilot_passed_for_this_state'] is not True or review['material_admission']):
            raise ValueError('Native state must have the exact independent/resource qualification and scope')
        if sha256_file(review['namespace'])!=review['namespace_sha256']:raise ValueError('State namespace changed')
        reader.verify_all(review['source_sha256']);reader.recheck()
        raw=json.loads((Path(spec['case'])/'result.json').read_text())
        native=json.loads(reader.path(raw['native_record_path']).read_text())
        energy=json.loads(reader.path(raw['active_record_path']).read_text())
        if native['identity']['protocol']!=NATIVE_THERMAL or energy['protocol']!=ACTIVE_SP:
            raise ValueError('New molecular states do not share the declared active/native protocols')
        if any(abs(v-raw['G_by_cutoff_eV'][c])>1e-8 for c,v in review['G_by_cutoff_eV'].items()):
            raise ValueError('Published state free energies differ from the reviewed backing calculation')
        partitions=[p for p in review['source_sha256'] if Path(p).parent.name=='projected_thermal' and p.endswith('.json')]
        if len(partitions)!=1:raise ValueError('One accepted projected partition per molecular state required')
        part=json.loads(Path(partitions[0]).read_text())
        numeric.append(dict(role='candidate',state=review['state'],charge=review['charge'],mult=review['mult'],
            product_inchikey=review['product_inchikey'],active_E_smd_eV=energy['e_smd_eV'],
            nuclear_by_cutoff_eV=part['nuclear_by_cutoff_eV'],spin_partition_eV=part['g_spin_partition_eV'],
            G_by_cutoff_eV=review['G_by_cutoff_eV'],active_record_sha256=sha256_file(reader.path(raw['active_record_path'])),
            partition_sha256=sha256_file(partitions[0]),native_protocol=native['identity']['protocol'],active_protocol=energy['protocol']))
        reviews.append(review)
    reference=json.loads(reference_path.read_text());reader=ArchivedSources(ROOT,archive,relative_to_prior=True)
    if reference['status']!='TEMPO_local_oxidation_native_partition_independently_reconstructed':
        raise ValueError('The selected matched native Fc source has not been independently reconstructed')
    reader.verify_all(reference['source_sha256']);reader.recheck()
    fc={name:reference['states']['ferrocene__'+name] for name in ('neu','ox')}
    definition=importlib.util.spec_from_file_location('independent_relocated_Fc_partition',archive/'scripts/validation/audit/native_partition.py')
    audit=importlib.util.module_from_spec(definition);definition.loader.exec_module(audit)
    for name in ('neu','ox'):
        base=str(ROOT/'calcs/dft/ferrocene'/name)
        candidates=[p for p in reference['source_sha256'] if p.startswith(base+'/native_thermal/') and Path(p).parent.name=='native_thermal' and p.endswith('.json')]
        if len(candidates)!=1 or json.loads(reader.path(candidates[0]).read_text())['identity']['protocol']!=NATIVE_THERMAL:
            raise ValueError('Fc and candidate native protocols differ')
        energies=[p for p in reference['source_sha256'] if p.startswith(base+'/sp/') and Path(p).parent.name=='sp' and p.endswith('.json')]
        if len(energies)!=1:raise ValueError('One checked active Fc energy required')
        energy=json.loads(reader.path(energies[0]).read_text())
        partition=audit.read_partition(reader.path(candidates[0]),reader.path(base+'/opt.xyz'))
        if energy['protocol']!=ACTIVE_SP or partition is None:raise ValueError('Matched Fc energy/nuclear reconstruction failed')
        reconstructed={str(c):energy['e_smd_eV']+partition['nuclear_by_cutoff_eV'][str(c)]
            +partition['g_thermal_total_eV']-partition['nuclear_by_cutoff_eV']['50'] for c in (25,50,100)}
        if any(abs(value-fc[name]['G_by_cutoff_eV'][c])>1e-8 for c,value in reconstructed.items()):
            raise ValueError('Fc source aggregate differs from its independent raw partition reconstruction')
        numeric.append(dict(role='reference',state=name,charge=fc[name]['charge'],mult=fc[name]['multiplicity'],
            active_E_smd_eV=energy['e_smd_eV'],nuclear_by_cutoff_eV=partition['nuclear_by_cutoff_eV'],
            spin_partition_eV=partition['g_thermal_total_eV']-partition['nuclear_by_cutoff_eV']['50'],
            G_by_cutoff_eV=reconstructed,native_protocol=NATIVE_THERMAL,active_protocol=energy['protocol']))
    values=ladder(reviews,fc)
    if len({r['candidate_alias'] for r in reviews})!=1:raise ValueError('Native states belong to different candidate aliases')
    pool=json.loads(registry_path.read_text())['candidates'];matches=[r for r in pool if r['inchikey']==reviews[0]['product_inchikey']]
    if len(matches)!=1 or matches[0]['route_status']!='proceed':raise ValueError('One exact currently route-admitted registry graph required')
    candidate=matches[0];policy=json.loads(policy_path.read_text())
    axes={axis['id']:dict(product_key=candidate['inchikey'],scope=axis['scope'],status='missing',reason='Exact application criterion/evidence remains unresolved') for axis in policy['required_axes']}
    source=[dict(path=str(status_path),sha256=sha256_file(status_path))]
    axes['attachment_route'].update(status='qualified',value='proceed',sources=[dict(path=str(registry_path),sha256=sha256_file(registry_path))])
    axes['redox_window'].update(status='partial',computed_events=values['events'],sources=source,
        reason='Checked capped-model potentials; exact-graft prediction interval and application cutoff unresolved')
    axes['state_stability'].update(status='partial',sources=source,
        reason='Declared-spin monomer calculations and disproportionation; global spin ordering and degradation/dimerization unresolved')
    if electrode_pool not in ('anolyte','catholyte'):raise ValueError('Declare the assessed electrode pool explicitly')
    decision=assess_candidate(dict(id=reviews[0]['candidate_alias'],product_key=candidate['inchikey'],role='candidate',electrode_pools=[electrode_pool],axes=axes),policy)
    result=dict(status='exact_molecular_native_redox_published_with_field_QC',candidate_alias=reviews[0]['candidate_alias'],
        product_smiles=candidate['product_smiles'],product_inchikey=candidate['inchikey'],registry_graph_key=candidate['graph_key'],
        model_scope='capped local graft; no explicit polymer or counterions',redox=values,
        grafting_free_energy_kJ_mol=None,occupancy_nonadditivity_kJ_mol=None,electronic_coupling_eV=None,
        directional_solution_relaxation_eV=None,all_axis_decision=decision,material_admission=False,
        global_conformer_convergence_proven=False,experimental_exact_graft_accuracy_validated=False,
        native_state_reviews=[dict(path=r['review'],sha256=r['review_sha256']) for r in status['entries']],
        source_bindings={str(p):sha256_file(p) for p in (status_path,reference_path,registry_path,policy_path)},
        reference_scope='Matched native Fc; spin-only partition; no empirical offset or orbital-degeneracy correction')
    output.mkdir(parents=True);inputs=output/'inputs';inputs.mkdir()
    for path,name in [(status_path,'state_collection.json'),(reference_path,'native_reference.json'),(policy_path,'all_axes_policy.json')]:shutil.copyfile(path,inputs/name)
    for review,row in zip(reviews,status['entries']):shutil.copyfile(row['review'],inputs/(review['state']+'_review.json'))
    atomic_write_json(output/'numeric_inputs.json',dict(status='exact_candidate_and_matched_reference_native_components',states=numeric))
    primary=[event['potential_V_vs_Fc_by_cutoff']['50'] for event in values['events']]
    names=['First','Second','Third','Fourth']
    if len(primary)>len(names):raise ValueError('Explicit TeX macro naming required for a longer redox ladder')
    macros={'Molecular'+names[i]+'Reduction':f'{value:+.3f}' for i,value in enumerate(primary)}
    if 'spacing_V_by_cutoff' in values:
        macros.update(MolecularSpacing=f"{values['spacing_V_by_cutoff']['50']:.3f}",
            MolecularDisproportionation=f"{values['disproportionation_kJ_mol_by_cutoff']['50']:.2f}")
    (output/'values.tex').write_text(''.join('\\newcommand{\\'+key+'}{'+value+'}\n' for key,value in macros.items()))
    atomic_write_json(output/'candidate.json',result);return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('states','reference','archive','registry','policy','out'):p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--electrode-pool',choices=('anolyte','catholyte'),required=True)
    a=p.parse_args()
    from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity
    import os
    with allocation_guard():host_capacity(1,minimum_ram_gib=4);os.sched_setaffinity(0,cpu_allocation(1))
    result=publish(a.states.resolve(),a.reference.resolve(),a.archive.resolve(),a.registry.resolve(),a.policy.resolve(),a.out.resolve(),a.electrode_pool)
    print(json.dumps(dict(status=result['status'],events=result['redox']['events'],all_axis_decision=result['all_axis_decision']['status'])))
