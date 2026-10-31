#!/usr/bin/env python3
"""Join a checked solution descriptor to an exact candidate without filling gaps."""
import argparse
import json
import math
from pathlib import Path
import shutil
import time

from rdkit import Chem

from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json
from redox.screening.feasibility import assess_candidate


def assemble(review,native,registry,policy,bindings):
    if (review['status']!='exact_TEMPO_directional_solution_descriptor_independently_checked'
            or not review['summary']['calculation_verified'] or not review['summary']['endpoint_model_verified']
            or review['candidate_id']!='tempo' or native['product_inchikey']!=review['product_inchikey']
            or Chem.MolToSmiles(Chem.MolFromSmiles(native['product_smiles']))!=review['product_smiles']):
        raise ValueError('Only the exact independently checked local candidate may be published')
    matches=[r for r in registry['candidates'] if r['inchikey']==review['product_inchikey']]
    if len(matches)!=1 or matches[0]['product_smiles']!=review['product_smiles']:
        raise ValueError('Descriptor and unified-registry graphs differ')
    if matches[0]['route_status']!='proceed':raise ValueError('This admitted-case publication requires the current exact route pass')
    if bindings['native']['sha256'] not in review['source_sha256'].values():
        raise ValueError('The published redox source must be the one bound by this pair review')
    record=matches[0];product_key=record['inchikey'];summary=review['summary']
    reduction,oxidation=summary['lambda_reduction_eV'],summary['lambda_oxidation_eV']
    if not all(math.isfinite(v) and v>0 for v in (reduction,oxidation)):
        raise ValueError('This checked positive-relaxation publication requires positive finite directions')
    relative=abs(reduction-oxidation)/((reduction+oxidation)/2)
    if (abs(relative-summary['directional_relative_asymmetry'])>1e-12
            or abs((reduction+oxidation)/2-summary['symmetric_single_body_approximation_eV'])>1e-12
            or any(abs(r['total_solution_relaxation_eV']-v)>1e-8 for r,v in zip(review['directional'],(reduction,oxidation)))):
        raise ValueError('Directional or summary arithmetic differs')
    sources=lambda name:[dict(path=bindings[name]['path'],sha256=bindings[name]['sha256'])]
    evidence={axis['id']:dict(product_key=product_key,scope=axis['scope'],status='missing',
                             reason='Required exact-case criterion/evidence remains unresolved') for axis in policy['required_axes']}
    evidence['attachment_route'].update(status='qualified',value='proceed',sources=sources('review'),
        context=dict(route='4-hydroxy-TEMPO attachment to preformed Merrifield resin',base='NaH',solvent='DMF or THF'))
    evidence['redox_window'].update(status='partial',reason='Checked capped-model potential; candidate-specific prediction interval and material transfer unresolved',
        model_value_V_vs_Fc=native['potential_V_vs_Fc_by_rotor_cutoff']['50'],sources=sources('native'))
    evidence['charge_transfer'].update(status='partial',reason='Directional composite solution descriptors checked; application cutoff, ensemble and rate model unresolved',
        directional_eV=dict(reduction=summary['lambda_reduction_eV'],oxidation=summary['lambda_oxidation_eV']),
        sources=sources('review'))
    candidate=dict(id='tempo',product_key=product_key,role='candidate',electrode_pools=['catholyte'],axes=evidence)
    decision=assess_candidate(candidate,policy)
    symmetric_ok=relative<=summary['equal_curvature_diagnostic_tolerance']
    if symmetric_ok!=('equal_curvature_approximation_unverified' not in summary['physical_flags']):
        raise ValueError('Symmetric-model flag contradicts its declared directional diagnostic')
    return dict(status='exact_candidate_solution_descriptor_published_with_field_QC',candidate_id='tempo',
        product_smiles=review['product_smiles'],product_inchikey=product_key,registry_graph_key=record['graph_key'],
        registry_aliases=record['aliases'],model_scope='capped local graft; no explicit polymer or counterions',
        reaction_context='NaH/DMF or THF attachment precedent; no reaction free energy computed',
        operating_model='native SMD MeCN frozen-slow electronic relaxation at geometry-level minima',
        redox_potential_V_vs_Fc=native['potential_V_vs_Fc_by_rotor_cutoff']['50'],
        solution_relaxation_eV=dict(reduction=summary['lambda_reduction_eV'],oxidation=summary['lambda_oxidation_eV']),
        components=review['directional'],directional_relative_asymmetry=summary['directional_relative_asymmetry'],
        equal_curvature_diagnostic_tolerance=summary['equal_curvature_diagnostic_tolerance'],
        qualified_symmetric_single_body_eV=summary['symmetric_single_body_approximation_eV'] if symmetric_ok else None,
        diagnostic_average_eV=summary['symmetric_single_body_approximation_eV'],
        grafting_free_energy_kJ_mol=None,occupancy_nonadditivity_kJ_mol=None,
        unwanted_association_free_energy_kJ_mol=None,electronic_coupling_eV=None,
        field_QC=dict(redox='independently reconstructed native local-model potential',
                      directional_solution_relaxation='independently checked composite electronic descriptor',
                      symmetric_solution_relaxation='qualified_diagnostic' if symmetric_ok else 'equal_curvature_diagnostic_failed',
                      grafting='not_computed',occupancy='not_computed',coupling='not_computed'),
        all_axis_decision=decision,experimental_rate_inferred=False,material_admission=False,
        source_bindings=bindings,created=time.time())


def publish(review_path,native_path,registry_path,policy_path,output):
    if output.exists():raise FileExistsError('Preserve earlier candidate publications')
    paths=dict(review=review_path,native=native_path,registry=registry_path,policy=policy_path)
    bindings={k:dict(path=str(p.resolve()),sha256=sha256_file(p)) for k,p in paths.items()}
    data={k:json.loads(p.read_text()) for k,p in paths.items()}
    for path,digest in data['review']['source_sha256'].items():
        if sha256_file(path)!=digest:raise ValueError('Checked solution backing record changed')
    raw_paths=[Path(p) for p in data['review']['source_sha256'] if Path(p).name=='result.json' and Path(p).parent.name=='pair_results']
    if len(raw_paths)!=1:raise ValueError('One bound raw solution-pair result required')
    raw=json.loads(raw_paths[0].read_text())
    if raw['directional']!=data['review']['directional'] or raw['identity']['method']!=data['review']['method']:
        raise ValueError('Published directional results differ from their bound raw calculation')
    pair=raw['identity']['pair'];geometry_names={raw['identity']['input_sha256'][s['geometry']]:s['name'] for s in pair['states']}
    numeric=[]
    for path in sorted(raw_paths[0].parent.glob('scf_*.json')):
        digest=data['review']['source_sha256'].get(str(path))
        if digest is None or sha256_file(path)!=digest:raise ValueError('Unbound SCF in the published pair')
        cell=json.loads(path.read_text());key=cell['key'];state=key['state']
        numeric.append(dict(state=state['name'],charge=state['charge'],mult=state['mult'],
            at=geometry_names[key['geometry_sha256']],geometry_sha256=key['geometry_sha256'],
            reference=None if key['reference'] is None else key['reference']['name'],
            reference_density_sha256=key['reference_density_sha256'],e_Ha=cell['e_Ha'],e_eV=cell['e_eV'],
            converged=cell['converged'],scf_guesses=cell['scf_guesses'],s_squared=cell['s2'],
            electron_trace=cell['electron_trace'],density_matrix_sha256=cell['density_matrix_sha256'],
            density_artifact_sha256=cell['density_artifact_sha256'],record_sha256=digest))
    if len(numeric)!=6:raise ValueError('Exactly six checked energy/density cells required')
    output.mkdir(parents=True);snapshot=output/'inputs';snapshot.mkdir()
    for key,path in paths.items():shutil.copyfile(path,snapshot/(key+'.json'))
    result=assemble(data['review'],data['native'],data['registry'],data['policy'],bindings)
    if any(sha256_file(b['path'])!=b['sha256'] for b in bindings.values()):raise ValueError('Candidate publication inputs changed')
    atomic_write_json(output/'numeric_inputs.json',dict(status='six_checked_solution_pair_numeric_cells',
        cells=numeric,method=raw['identity']['method'],input_sha256=raw['identity']['input_sha256'],
        raw_result_sha256=sha256_file(raw_paths[0]),raw_result_path=str(raw_paths[0])))
    atomic_write_json(output/'candidate.json',result);return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ('review','native','registry','policy','out'):parser.add_argument('--'+name,type=Path,required=True)
    args=parser.parse_args();result=publish(args.review,args.native,args.registry,args.policy,args.out)
    print(json.dumps(dict(status=result['status'],candidate=result['candidate_id'],relaxation=result['solution_relaxation_eV'],
                         symmetric_scalar=result['qualified_symmetric_single_body_eV'],all_axis_decision=result['all_axis_decision']['status'])))
