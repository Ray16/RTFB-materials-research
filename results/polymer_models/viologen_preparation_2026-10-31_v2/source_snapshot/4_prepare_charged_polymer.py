#!/usr/bin/env python3
"""Prepare route-bound ionic occupancy graphs and audited initial coordinates.

This command never starts UMA, xTB, SCF, QM sampling or a GPU reservation.
"""
import argparse
import importlib.util
import json
import os
from pathlib import Path
import shutil

from rdkit import Chem,rdBase
from rdkit.Chem import AllChem,Descriptors,rdDetermineBonds

from redox.core.common import ROOT
from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json
from redox.build.geometry import _skeleton,geometry_stereo_status
from redox.build.ionic_polymer import charged_panel
from redox.build.oligomer_topology import chain_topology
from redox.build.polymer_contacts import pi_core_groups,contact_metrics
from redox.validation.protonation import hydrogen_topology_status


def geometry(case,fragment,settings,output,seed):
    mol=Chem.AddHs(Chem.MolFromSmiles(case['smiles']));parameters=AllChem.ETKDGv3()
    parameters.randomSeed=seed;parameters.enforceChirality=True;parameters.numThreads=1
    parameters.maxIterations=settings['embedding_max_iterations']
    parameters.trackFailures=True
    ids=list(AllChem.EmbedMultipleConfs(mol,numConfs=settings['embeds_per_graph'],params=parameters))
    attempts=[dict(profile='ETKDGv3',seed=seed,conformers=len(ids),failure_counts=list(parameters.GetFailureCounts()))]
    if not ids and settings['random_coordinate_fallback']:
        parameters.useRandomCoords=True
        ids=list(AllChem.EmbedMultipleConfs(mol,numConfs=settings['embeds_per_graph'],params=parameters))
        attempts.append(dict(profile='ETKDGv3_random_coordinates',seed=seed,conformers=len(ids),failure_counts=list(parameters.GetFailureCounts())))
    atomic_write_json(output/'embedding.json',dict(attempts=attempts,settings=settings,rdkit_version=rdBase.rdkitVersion))
    if not ids:raise ValueError('No initial embedding obtained after the declared bounded attempts')
    props=AllChem.MMFFGetMoleculeProperties(mol)
    if props is None:raise ValueError('No declared MMFF parameters for this graph')
    results=AllChem.MMFFOptimizeMoleculeConfs(mol,numThreads=1,maxIters=settings['MMFF_max_iterations'])
    topology=chain_topology(case['smiles'],case['repeat_units'],case['reactive_sites'])
    trials=[];selected=None
    for cid in sorted(ids,key=lambda i:results[i][1]):
        xyz=output/f'trial_{cid}.xyz';Chem.MolToXYZFile(mol,str(xyz),confId=cid)
        connectivity=hydrogen_topology_status(case['smiles'],xyz);stereo=geometry_stereo_status(case['smiles'],xyz)
        row=dict(conformer_id=cid,MMFF_status=results[cid][0],MMFF_energy_kcal_mol=results[cid][1],connectivity=connectivity,stereochemistry=stereo,xyz=str(xyz),xyz_sha256=sha256_file(xyz))
        trials.append(row)
        if selected is None and connectivity=='ok' and stereo in ('ok','not_applicable'):selected=row
    if selected is None:raise ValueError('Every initial geometry failed connectivity/stereo review')
    observed=Chem.MolFromXYZFile(selected['xyz']);rdDetermineBonds.DetermineConnectivity(observed,charge=0)
    expected=Chem.AddHs(Chem.MolFromSmiles(case['smiles']));mapping=_skeleton(observed).GetSubstructMatch(_skeleton(expected))
    if len(mapping)!=expected.GetNumAtoms():raise ValueError('Full graph cannot map onto initial XYZ')
    backbone={mapping[i] for i in topology['backbone_atoms']};degrees=[sum(a.GetIdx() in backbone for a in observed.GetAtomWithIdx(i).GetNeighbors()) for i in backbone]
    edges=sum(degrees)//2
    if len(backbone)!=2*case['repeat_units'] or edges!=2*case['repeat_units']-1 or degrees.count(1)!=2 or max(degrees)>2:
        raise ValueError('Actual XYZ backbone is not the intended open path')
    molecule=Chem.MolFromSmiles(case['smiles']);groups=pi_core_groups(molecule,fragment)
    if len(groups)!=len(case['occupied_sites']):raise ValueError('Motif ring-core count differs from occupied sites')
    positions=observed.GetConformer().GetPositions();xyz_groups=[[mapping[i] for i in group] for group in groups]
    # Distances/planes describe only this initial embedding, not its equilibrium
    # structure, interaction strength, coupling or conformer population.
    contacts=contact_metrics(positions,xyz_groups)
    graph_mass=Descriptors.MolWt(molecule)
    for bound in case['composition_bounds'].values():
        if abs(bound['polymer_cation_mass_g_mol']-graph_mass)>1e-8:raise ValueError('Analytic charged-polymer mass differs from its actual graph')
    return dict(status='initial_connectivity_stereo_and_backbone_checked',selected=selected,trials=trials,embedding_attempts=attempts,
        topology=topology,XYZ_backbone_atoms=sorted(backbone),XYZ_backbone_carbon_count=len(backbone),
        XYZ_backbone_CC_bonds=edges,XYZ_chain_ends=2,XYZ_core_groups=xyz_groups,initial_contact_metrics=contacts,
        solvent_optimized=False,conformer_sampling_qualified=False,interaction_strength_computed=False,native_G_qualified=False)


def prepare(config_path,registry_path,routes_path,assessment_path,output):
    if output.exists():raise FileExistsError('Preserve earlier preparation outputs')
    config=json.loads(config_path.read_text());registry=json.loads(registry_path.read_text())['candidates']
    if config['quantum_calculation_admitted'] is not False or config.get('occupancy_weights') is not None or config['experimental_loading_known'] is not False:
        raise ValueError('This command prepares explicitly conditional, unqualified charged models only')
    if len(config['reactive_sites'])!=len(set(config['reactive_sites'])):raise ValueError('Reactive sites must be unique')
    matches=[r for r in registry if config['candidate_alias'] in r['aliases']]
    if len(matches)!=1:raise ValueError('One exact registered local graft required')
    candidate=matches[0];declared=next(d for d in candidate['declarations'] if d['id']==config['candidate_alias']);fragment=declared['fragment']
    definition=importlib.util.spec_from_file_location('charged_preparation_route_gate',ROOT/'scripts/pipeline/screen_grafting_routes.py');gate=importlib.util.module_from_spec(definition);definition.loader.exec_module(gate)
    route=gate.screen(dict(id=config['candidate_alias'],frag=fragment,product_smiles=candidate['product_smiles']),json.loads(routes_path.read_text())['candidates'],json.loads(assessment_path.read_text()))
    if not route['proceed_to_full_redox'] or not route.get('exact_route_binding_verified'):raise ValueError('The exact local attachment route must pass prospectively')
    all_cases=[];all_cycles=[]
    for name,stereo in config['configurations'].items():
        spec=dict(config,backbone_stereo_tags=stereo);cases,cycles,ionic=charged_panel(spec,fragment)
        for case in cases:case.update(configuration=name,key=name+'__'+case['key'])
        for cycle in cycles:
            cycle['coefficients']={name+'__'+k:v for k,v in cycle['coefficients'].items()};cycle['configuration']=name
        all_cases.extend(cases);all_cycles.extend(cycles)
    inputs=[config_path,registry_path,routes_path,assessment_path,Path(__file__),ROOT/'src/redox/build/polymer_cases.py',ROOT/'src/redox/build/ionic_polymer.py',ROOT/'src/redox/build/oligomer_topology.py',ROOT/'src/redox/build/geometry.py',ROOT/'src/redox/build/polymer_contacts.py',ROOT/'src/redox/validation/protonation.py']
    hashes={str(p.resolve()):sha256_file(p) for p in inputs}
    output.mkdir(parents=True);reports=[];failures=[]
    snapshot=output/'source_snapshot';snapshot.mkdir()
    copies={}
    for index,path in enumerate(inputs):
        target=snapshot/(str(index)+'_'+path.name);shutil.copyfile(path,target)
        if sha256_file(target)!=hashes[str(path.resolve())]:raise ValueError('Source changed during preparation snapshot')
        copies[str(path.resolve())]=str(target.relative_to(output))
    atomic_write_json(output/'plan.json',dict(config=config,candidate_alias=config['candidate_alias'],route_gate=route,cases=all_cases,cycles=all_cycles,source_sha256=hashes,snapshot_files=copies,rdkit_version=rdBase.rdkitVersion,quantum_jobs_submitted=0))
    for index,case in enumerate(all_cases):
        folder=output/case['key'];folder.mkdir()
        atomic_write_json(folder/'case.json',case)
        atomic_write_json(output/'progress.json',dict(prepared=len(reports),failed=failures,active_case=case['key'],expected=len(all_cases),quantum_jobs_submitted=0))
        try:
            audit=geometry(case,fragment,config['initial_geometry'],folder,config['initial_geometry']['seed']+index)
        except (ValueError,RuntimeError) as error:
            failure=dict(case=case['key'],status='initial_geometry_failed',error_type=type(error).__name__,error=str(error),quantum_jobs_submitted=0)
            atomic_write_json(folder/'failure.json',failure);failures.append(failure)
        else:
            atomic_write_json(folder/'geometry_review.json',audit)
            reports.append(dict(case=case,geometry_review=str(folder/'geometry_review.json'),geometry_review_sha256=sha256_file(folder/'geometry_review.json'),initial_geometry=audit['selected']['xyz'],initial_geometry_sha256=audit['selected']['xyz_sha256']))
        atomic_write_json(output/'progress.json',dict(prepared=len(reports),failed=failures,active_case=None,expected=len(all_cases),quantum_jobs_submitted=0))
    if any(sha256_file(Path(path))!=value for path,value in hashes.items()):raise ValueError('Preparation source changed during geometry construction')
    result=dict(status='initial_geometry_partial' if failures else 'exact_route_bound_ionic_polymer_graphs_and_initial_geometries_prepared',candidate_alias=config['candidate_alias'],local_graft_inchikey=candidate['inchikey'],local_graft_smiles=candidate['product_smiles'],fragment=fragment,config=config,route_gate=route,ionic_definition=ionic,cases=reports,cycles=all_cycles,failures=failures,expected_graphs=len(all_cases),
        quantum_jobs_submitted=0,quantum_calculation_admitted=False,material_admission=False,
        charged_polymer_solver_and_resource_qualification_complete=False,ion_association_or_exchange_computed=False,
        native_G_eV=None,experimental_loading_inferred=False,source_sha256=hashes,snapshot_files=copies,rdkit_version=rdBase.rdkitVersion,
        files_sha256={str(p.relative_to(output)):sha256_file(p) for p in output.rglob('*') if p.is_file()})
    atomic_write_json(output/'preparation.json',result);return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('config','registry','routes','assessment','out'):p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args()
    from redox.workflows.resources import allocation_guard,host_capacity,cpu_allocation
    with allocation_guard():host_capacity(1,minimum_ram_gib=4);os.sched_setaffinity(0,cpu_allocation(1))
    result=prepare(a.config.resolve(),a.registry.resolve(),a.routes.resolve(),a.assessment.resolve(),a.out.resolve());print(json.dumps(dict(status=result['status'],graphs=len(result['cases']),cycles=len(result['cycles']),quantum_jobs_submitted=0)))
    if result['failures']:raise SystemExit(1)
