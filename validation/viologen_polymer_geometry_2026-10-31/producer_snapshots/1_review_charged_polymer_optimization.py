#!/usr/bin/env python3
"""Independently bind a converged ionic geometry to its final SCF and graph.

This checks stationarity and retained connectivity/stereochemistry. A Hessian,
conformer ensemble, interaction energy and material eligibility remain separate.
"""
import argparse
import importlib.util
import json
import math
import os
from pathlib import Path

import h5py
import numpy as np
from ase.io import read as read_atoms
from pyscf import gto
from pyscf.lib import param
from rdkit import Chem
from rdkit.Chem import rdDetermineBonds

from redox.build.geometry import _skeleton,geometry_stereo_status
from redox.build.oligomer_topology import chain_topology
from redox.core.common import ROOT
from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json
from redox.workflows.molecular_native_namespace import verify


def physical_checks(plan,returned,final,xyz_coordinates,checkpoint_coordinates,checkpoint_energy):
    n=plan['atoms'];coordinates=np.asarray(final['coordinates_Bohr'],dtype=float)
    force=np.asarray(final['gradient_Ha_Bohr'],dtype=float)
    if (coordinates.shape!=(n,3) or force.shape!=(n,3) or not np.isfinite(coordinates).all()
            or not np.isfinite(force).all() or not math.isfinite(final['E_Ha'])
            or not returned['converged'] or not final['optimizer_converged'] or not final['SCF_converged']
            or final['charge']!=plan['charge'] or final['mult']!=plan['mult']
            or returned['charge']!=plan['charge'] or returned['mult']!=plan['mult']
            or abs(final['spin'].get('spin_contam',math.inf))>.1 or final['native_G_qualified']):
        raise ValueError('Finite, converged returned-geometry SCF/force/state records required')
    for actual in (xyz_coordinates,np.asarray(returned['coordinates_A']),checkpoint_coordinates):
        actual=np.asarray(actual)
        if actual.shape!=(n,3) or not np.allclose(coordinates*param.BOHR,actual,rtol=0,atol=1e-8):
            raise ValueError('Optimizer, final evaluation, XYZ and actual SCF checkpoint geometries differ')
    if not math.isfinite(checkpoint_energy) or abs(final['E_Ha']-checkpoint_energy)>1e-10:raise ValueError('Returned geometry has a different checkpoint energy')
    atom_norms=np.linalg.norm(force,axis=1);maximum=float(atom_norms.max());rms=float(np.sqrt(np.mean(atom_norms**2)))
    if (maximum>=plan['optimizer_settings']['convergence_gmax'] or rms>=plan['optimizer_settings']['convergence_grms']
            or not math.isclose(float(abs(force).max()),final['gradient_max_Ha_Bohr'],abs_tol=1e-14,rel_tol=0)):
        raise ValueError('Fresh final gradient fails the declared per-atom norm stationarity limits')
    return dict(maximum_atom_gradient_Ha_Bohr=maximum,RMS_atom_gradient_Ha_Bohr=rms,
                maximum_Cartesian_gradient_Ha_Bohr=float(abs(force).max()))


def graph_checks(smiles,initial_xyz,final_xyz):
    expected=Chem.AddHs(Chem.MolFromSmiles(smiles));observed=[]
    for path in (initial_xyz,final_xyz):
        mol=Chem.MolFromXYZFile(str(path))
        if mol is None:raise ValueError('Unreadable charged-polymer geometry')
        rdDetermineBonds.DetermineConnectivity(mol,charge=0);observed.append(mol)
    for mol in observed:
        if (mol.GetNumAtoms()!=expected.GetNumAtoms() or mol.GetNumBonds()!=expected.GetNumBonds()
                or len(_skeleton(mol).GetSubstructMatch(_skeleton(expected)))!=expected.GetNumAtoms()):
            raise ValueError('Full inferred covalent graph differs from the exact declared structure')
    bonds=lambda mol:{tuple(sorted((b.GetBeginAtomIdx(),b.GetEndAtomIdx()))) for b in mol.GetBonds()}
    if ([a.GetSymbol() for a in observed[0].GetAtoms()]!=[a.GetSymbol() for a in observed[1].GetAtoms()]
            or bonds(observed[0])!=bonds(observed[1])):
        raise ValueError('Indexed atom composition or connectivity changed during optimization')
    if geometry_stereo_status(smiles,final_xyz)!='ok':raise ValueError('Declared backbone stereochemistry did not survive')
    positions=observed[1].GetConformer().GetPositions();table=Chem.GetPeriodicTable();ratios=[]
    for bond in observed[1].GetBonds():
        i,j=bond.GetBeginAtomIdx(),bond.GetEndAtomIdx()
        radius=sum(table.GetRcovalent(observed[1].GetAtomWithIdx(k).GetAtomicNum()) for k in (i,j))
        ratios.append(float(np.linalg.norm(positions[i]-positions[j]))/radius)
    if min(ratios)<.75 or max(ratios)>1.30:raise ValueError('Implausible covalent bond length in the optimized geometry')
    return dict(connectivity_retained=True,stereochemistry_retained=True,
                minimum_bond_covalent_radius_ratio=min(ratios),maximum_bond_covalent_radius_ratio=max(ratios))


def review(namespace):
    spec=verify(namespace);case=Path(spec['case']);sources=Path(spec['sources']);used=[namespace]
    def read(path):used.append(path);return json.loads(path.read_text())
    plan=read(case/'plan.json');stage=read(case/'status.json')
    definition=importlib.util.spec_from_file_location('charged_geometry_closed_owners',ROOT/'scripts/pipeline/observe_archived_native_namespace.py')
    owners=importlib.util.module_from_spec(definition);definition.loader.exec_module(owners)
    placements=[p.parent for p in case.parent.glob('placement*/launch.json')]
    if len(placements)!=1:raise ValueError('One exact placement required')
    placement=placements[0];launch=read(placement/'launch.json');completion=read(placement/'completion.json')
    host=read(case/'host_guard/status.json');runtime=read(case/'tensor_runtime/status.json')
    owner=read(case/'residency_guard/owner.json');residency=read(case/'residency_guard/state.json')
    if (not plan.get('charged_polymer_optimization') or plan['role']!='charged_polymer_solution_geometry_pilot' or plan['native_G_qualified']
            or stage['status']!='charged_polymer_solution_geometry_finished_independent_review_pending'
            or owners.alive(stage['owner']) is not False or stage['native_G_qualified'] or stage['material_admission']
            or completion['returncode']!=0 or completion['budget_stop'] is not None
            or not owner['guard_stopped'] or not residency['own_residency_observed'] or residency['foreign_resident_PIDs']
            or host['status']!='effective_memory_guard_finished_profile_review_pending'
            or runtime['status']!='task_returned_under_checked_cutensor_runtime'
            or host['maximum_RSS_GiB']!=40 or host['free_RAM_bytes']<32*2**30
            or max(host['peak_tree_RSS_bytes'],host['self_VmHWM_bytes'])>40*2**30
            or not 0<stage['observed_peak_GPU_MiB']<=31744
            or not stage['owner']['pid']==owner['pid']==runtime['pid']==host['pid']
            or not stage['owner']['host']==owner['host']==host['host']==launch['host']
            or not stage['owner']['process_start']==owner['process_start']==host['process_start']
            or not owner['CPU_affinity']==runtime['CPU_affinity']==host['CPU_affinity']==launch['CPU_affinity']):
        raise ValueError('Normal closed optimization, exclusive ownership and bounded resources required')
    route=read(sources/'prospective_route.json')
    if not route['proceed_to_full_redox'] or not route['exact_route_binding_verified']:raise ValueError('Prospective exact route must pass')
    returned=read(case/'optimizer_return.json');final=read(case/'final_evaluation.json')
    xyz=case/'returned.xyz';initial=sources/'input.xyz';checkpoint=case/'final_SCF.h5';used.extend([xyz,initial,checkpoint])
    atoms=read_atoms(xyz)
    with h5py.File(checkpoint,'r') as handle:
        mol=gto.loads(handle['mol'][()]);energy=float(handle['scf/e_tot'][()]);occupations=handle['scf/mo_occ'][()]
    if (mol.charge!=plan['charge'] or mol.spin+1!=plan['mult'] or mol.basis!='def2-svp'
            or not np.isfinite(occupations).all() or (occupations<0).any() or (occupations>2).any()
            or abs(float(occupations.sum())-mol.nelectron)>1e-10
            or returned['elements']!=atoms.get_chemical_symbols()
            or [mol.atom_symbol(i) for i in range(mol.natm)]!=atoms.get_chemical_symbols()):
        raise ValueError('Actual final SCF composition, charge, occupations or basis differ')
    physical=physical_checks(plan,returned,final,atoms.positions,mol.atom_coords(unit='Angstrom'),energy)
    geometry=graph_checks(plan['prepared_case']['smiles'],initial,xyz)
    evaluations=sorted((case/'geometry_evaluations').glob('*.json'))
    if len(evaluations)!=returned['geometry_evaluations'] or len(evaluations)!=stage['geometry_evaluations']:
        raise ValueError('Incomplete durable geometry evaluation history')
    for index,path in enumerate(evaluations,1):
        row=read(path)
        if (path.stem!=f'{index:04d}' or not row['SCF_converged'] or row['charge']!=plan['charge'] or row['mult']!=plan['mult']
                or row['elements']!=atoms.get_chemical_symbols() or row['owner']!=stage['owner']
                or np.asarray(row['coordinates_Bohr']).shape!=(plan['atoms'],3)
                or np.asarray(row['gradient_Ha_Bohr']).shape!=(plan['atoms'],3)
                or not np.isfinite(row['coordinates_Bohr']).all() or not np.isfinite(row['gradient_Ha_Bohr']).all()
                or not math.isfinite(row['E_Ha'])):raise ValueError('Invalid saved geometry evaluation')
    preparation=read(sources/'preparation.json');topology=chain_topology(plan['prepared_case']['smiles'],preparation['config']['repeat_units'],preparation['config']['reactive_sites'])
    return dict(status='charged_polymer_solution_stationary_geometry_independently_checked',namespace=str(namespace),
        graph_key=plan['prepared_case']['key'],candidate_alias=plan['exact_candidate_alias'],charge=plan['charge'],mult=plan['mult'],atoms=plan['atoms'],
        geometry_sha256=sha256_file(xyz),E_Ha=final['E_Ha'],geometry_evaluations=len(evaluations),**physical,**geometry,
        backbone_carbons=topology['backbone_carbon_count'],backbone_CC_bonds=topology['backbone_CC_bonds'],chain_ends=topology['chain_ends'],
        peak_GPU_MiB=stage['observed_peak_GPU_MiB'],solvent_model=plan['solvent_model'],
        Hessian_minimum_check_complete=False,conformer_sampling_qualified=False,interaction_strength_computed=False,
        native_G_qualified=False,material_admission=False,source_sha256={str(p):sha256_file(p) for p in used},producer_sha256=sha256_file(__file__))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--spec',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    from redox.workflows.resources import allocation_guard,host_capacity,cpu_allocation
    with allocation_guard():host_capacity(1,minimum_ram_gib=4);os.sched_setaffinity(0,cpu_allocation(1))
    if a.out.exists():raise FileExistsError('Preserve earlier optimization reviews')
    result=review(a.spec.resolve());atomic_write_json(a.out,result)
    print(json.dumps({k:result[k] for k in ('status','graph_key','E_Ha','maximum_atom_gradient_Ha_Bohr')}))
