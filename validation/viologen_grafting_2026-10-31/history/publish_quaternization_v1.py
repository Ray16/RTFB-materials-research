#!/usr/bin/env python3
"""Assemble an exact quaternization from independently checked native records."""
import argparse
import importlib.util
import json
from pathlib import Path

from redox.build.grafting import grafting_reaction
from redox.core.common import ROOT
from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json
from redox.core.protocol import ACTIVE_SP
from redox.qm.native_thermal import NATIVE_THERMAL
from redox.properties.grafting import quaternization_delta_g
from redox.validation.molecular_native import MolecularSources
from redox.workflows.molecular_native_namespace import verify


def load(path,name):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module


def reviewed_component(path,reagent):
    review=json.loads(path.read_text())
    expected='reaction_component_native_G_independently_checked' if reagent else 'exact_molecular_native_G_independently_checked'
    if review['status']!=expected or review['resource_pilot_passed_for_this_state'] is not True:
        raise ValueError('Independent native and resource admission required for every molecular component')
    spec_path=Path(review['namespace']);spec=verify(spec_path)
    if sha256_file(spec_path)!=review['namespace_sha256']:raise ValueError('Reviewed namespace changed')
    reader=MolecularSources(spec);reader.verify_all(review['source_sha256']);reader.recheck()
    case=Path(spec['case']);plan=json.loads((case/'plan.json').read_text());raw=json.loads((case/'result.json').read_text())
    energy_path=reader.path(raw['active_record_path']);native_path=reader.path(raw['native_record_path'])
    energy=json.loads(energy_path.read_text());native=json.loads(native_path.read_text())
    if energy['protocol']!=ACTIVE_SP or native['identity']['protocol']!=NATIVE_THERMAL:
        raise ValueError('Electronic or native component protocol differs')
    if (review['charge'],review['mult'])!=(energy['charge'],energy['mult']):raise ValueError('Component state differs')
    if any(abs(v-raw['G_by_cutoff_eV'][c])>1e-8 for c,v in review['G_by_cutoff_eV'].items()):
        raise ValueError('Reviewed and produced component free energies differ')
    partitions=[reader.path(p) for p in review['source_sha256'] if Path(p).parent.name=='projected_thermal' and Path(p).suffix=='.json']
    if len(partitions)!=1:raise ValueError('One exact projected nuclear partition required')
    partition=json.loads(partitions[0].read_text())
    if any(abs(v-energy['e_smd_eV']-partition['nuclear_by_cutoff_eV'][c]-partition['g_spin_partition_eV'])>1e-8
           for c,v in review['G_by_cutoff_eV'].items()):raise ValueError('Electronic/nuclear/spin components do not sum')
    component=dict(smiles=plan['smiles'] if reagent else plan['product_smiles'],charge=review['charge'],mult=review['mult'],
        G_by_cutoff_eV=review['G_by_cutoff_eV'],native_protocol=native['identity']['protocol'],active_protocol=energy['protocol'],
        active_E_smd_eV=energy['e_smd_eV'],conformer_index=plan.get('conformer_index'),
        nuclear_by_cutoff_eV=partition['nuclear_by_cutoff_eV'],spin_partition_eV=partition['g_spin_partition_eV'],
        minimum_internal_frequency_cm=min(partition['frequencies_cm']),maximum_gradient_Ha_Bohr=native['gradient_max'],
        maximum_sampled_column_error=max(c['relative_error'] for c in native['derivative_checks']),
        native_sha256=sha256_file(native_path),partition_sha256=sha256_file(partitions[0]),
        final_geometry_sha256=review['source_sha256'][str(reader.path(Path(plan['target_directory'])/'opt.xyz'))])
    return component,{**review['source_sha256'],str(path.resolve()):sha256_file(path),str(spec_path):sha256_file(spec_path)},plan


def publish(product_review,precursor_reviews,package,registry_path,recipe_path,output):
    if output.exists():raise FileExistsError('Preserve earlier reaction publications')
    recipe=json.loads(recipe_path.read_text());registry=json.loads(registry_path.read_text())['candidates']
    matches=[c for c in registry if recipe['route_candidate'] in c['aliases']]
    if len(matches)!=1 or matches[0]['route_status']!='proceed':raise ValueError('One exact route-admitted product required')
    candidate=matches[0];declaration=next(d for d in candidate['declarations'] if d['id']==recipe['route_candidate'])
    reaction=grafting_reaction(declaration['fragment'],candidate['product_smiles'])
    product,sources,product_plan=reviewed_component(product_review,False)
    if product_plan['candidate_alias']!=recipe['route_candidate']:raise ValueError('Product native state belongs to another design')
    precursors=[];indices=[]
    for path in precursor_reviews:
        component,binding,plan=reviewed_component(path,True)
        if (plan['source_id']!=recipe['source_id'] or plan['recipe']!=recipe
                or plan['route_candidate']!=recipe['route_candidate']):raise ValueError('Precursor recipe or route differs')
        indices.append(plan['conformer_index']);precursors.append(component);sources.update(binding)
    if sorted(indices)!=list(range(recipe['top_k'])):
        raise ValueError('All prospectively requested fixed precursor seeds must be checked before publishing')
    # The compact 14-state archive includes independently replayable local-site
    # and chloride records; it is never used to borrow a molecular correction.
    replay=load(package/'reproduce.py','quaternization_shared_reference_replay').replay(package,full_native=True)
    if not replay['nuclear_partition_reintegrated']:raise ValueError('Shared reference native replay required')
    manifest=json.loads((package/'manifest.json').read_text());shared={}
    for role,name in [('site','local'),('byproduct','chloride')]:
        row=next(r for r in manifest['states'] if r['id']==name)
        energy=json.loads((package/row['files']['energy']).read_text());part=json.loads((package/row['files']['partition']).read_text())
        if energy['protocol']!=ACTIVE_SP:raise ValueError('Shared reference energy protocol differs')
        if row['partition_kind']!='monatomic':
            native=json.loads((package/row['files']['native']).read_text())
            if native['identity']['protocol']!=NATIVE_THERMAL:raise ValueError('Shared reference native protocol differs')
        shared[role]=dict(smiles=row['smiles'],charge=row['charge'],mult=row['mult'],
            active_E_smd_eV=energy['e_smd_eV'],spin_partition_eV=0. if row['partition_kind']=='monatomic' else part['g_spin_partition_eV'],
            nuclear_by_cutoff_eV={str(c):part['G_thermal_eV'] for c in (25,50,100)} if row['partition_kind']=='monatomic' else part['nuclear_by_cutoff_eV'],
            G_by_cutoff_eV={str(c):energy['e_smd_eV']+(part['G_thermal_eV'] if row['partition_kind']=='monatomic'
                else part['nuclear_by_cutoff_eV'][str(c)]+part['g_spin_partition_eV']) for c in (25,50,100)})
    for relative,digest in manifest['files_sha256'].items():sources[str(package/relative)]=digest
    for path in (package/'manifest.json',package/'reproduce.py',package/'cached_reader.py',
                 package/'nuclear_equations.py',package/'projected_modes.py',registry_path,recipe_path,
                 Path(__file__).resolve(),ROOT/'src/redox/properties/grafting.py',ROOT/'src/redox/build/grafting.py'):
        sources[str(path.resolve())]=sha256_file(path)
    selected=min(range(len(precursors)),key=lambda i:precursors[i]['G_by_cutoff_eV']['50'])
    cycles=[quaternization_delta_g(reaction,dict(shared,product=product,precursor=p)) for p in precursors]
    result=dict(status='exact_local_quaternization_native_cycle_checked_finite_sampling',
        candidate_alias=recipe['route_candidate'],product_inchikey=candidate['inchikey'],reaction=reaction,
        selected_precursor_conformer_index=indices[selected],selection='lowest precursor G at 50 cm^-1 among three fixed seeds',
        **cycles[selected],components=dict(shared,product=product,precursor=precursors[selected]),
        precursor_conformers=precursors,per_precursor_cycle=cycles,
        finite_precursor_deltaG_span_kJ_mol=max(c['deltaG_kJ_mol_by_cutoff']['50'] for c in cycles)-min(c['deltaG_kJ_mol_by_cutoff']['50'] for c in cycles),
        source_sha256=sources,global_conformer_convergence_proven=False,
        experimental_grafting_accuracy_validated=False,shared_reference_replay=replay,
        scope='Finite capped-site cycle in MeCN at 298.15 K with separated chloride; no salt speciation, polymer penalty or reflux correction')
    if any(sha256_file(p)!=h for p,h in sources.items()):raise ValueError('A reaction source changed during assembly')
    output.mkdir(parents=True);atomic_write_json(output/'reaction.json',result)
    (output/'values.tex').write_text('\\newcommand{\\ViologenGraftDG}{'+f'{result["deltaG_kJ_mol_by_cutoff"]["50"]:.2f}'+'}\n')
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--product-review',type=Path,required=True);p.add_argument('--precursor-review',type=Path,action='append',required=True)
    for name in ('package','registry','recipe','out'):p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args()
    import os
    from redox.workflows.resources import allocation_guard,host_capacity,cpu_allocation
    with allocation_guard():host_capacity(1,minimum_ram_gib=4);os.sched_setaffinity(0,cpu_allocation(1))
    r=publish(a.product_review.resolve(),[f.resolve() for f in a.precursor_review],a.package.resolve(),a.registry.resolve(),a.recipe.resolve(),a.out.resolve())
    print(json.dumps(dict(status=r['status'],deltaG_kJ_mol_by_cutoff=r['deltaG_kJ_mol_by_cutoff'])))
