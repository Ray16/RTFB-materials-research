#!/usr/bin/env python3
"""Replay the ionic preparation's stored graphs, XYZ connectivity and masses.

This is an initial-structure audit, never a thermochemistry qualification.
"""
import argparse
from collections import Counter
import json
import math
import os
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors,rdDetermineBonds

from redox.build.geometry import _skeleton,geometry_stereo_status
from redox.build.oligomer_topology import chain_topology
from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json


def check(path):
    data=json.loads(path.read_text());root=path.parent
    if (data['status']!='exact_route_bound_ionic_polymer_graphs_and_initial_geometries_prepared'
            or data['failures'] or len(data['cases'])!=data['expected_graphs']
            or data['quantum_jobs_submitted'] or data['quantum_calculation_admitted']
            or data['material_admission'] or data['native_G_eV'] is not None):
        raise ValueError('A complete, explicitly unqualified preparation is required')
    for name,digest in data['files_sha256'].items():
        if sha256_file(root/name)!=digest:raise ValueError('Changed preparation file: '+name)
    for name,digest in data['source_sha256'].items():
        if sha256_file(root/data['snapshot_files'][name])!=digest:raise ValueError('Changed frozen source')
    config=data['config'];records=[];graphs={};table=Chem.GetPeriodicTable()
    for row in data['cases']:
        case=row['case'];key=case['key']
        if key in graphs:raise ValueError('Duplicate graph/configuration key')
        graph=Chem.AddHs(Chem.MolFromSmiles(case['smiles']));graphs[key]=graph
        topology=chain_topology(case['smiles'],config['repeat_units'],config['reactive_sites'])
        audit_path=Path(row['geometry_review']);xyz_path=Path(row['initial_geometry'])
        if sha256_file(audit_path)!=row['geometry_review_sha256'] or sha256_file(xyz_path)!=row['initial_geometry_sha256']:
            raise ValueError('Changed selected geometry or audit')
        audit=json.loads(audit_path.read_text())
        if audit['native_G_qualified'] or audit['solvent_optimized'] or audit['interaction_strength_computed']:
            raise ValueError('An initial geometry must not imply a qualified scientific result')
        observed=Chem.MolFromXYZFile(str(xyz_path));rdDetermineBonds.DetermineConnectivity(observed,charge=0)
        mapping=_skeleton(observed).GetSubstructMatch(_skeleton(graph))
        if len(mapping)!=graph.GetNumAtoms() or observed.GetNumAtoms()!=graph.GetNumAtoms() or observed.GetNumBonds()!=graph.GetNumBonds():
            raise ValueError('Full selected XYZ connectivity differs from the exact graph')
        if geometry_stereo_status(case['smiles'],xyz_path)!='ok':raise ValueError('Selected XYZ stereochemistry differs')
        positions=observed.GetConformer().GetPositions();ratios=[];attachments=[]
        for bond in graph.GetBonds():
            i,j=bond.GetBeginAtomIdx(),bond.GetEndAtomIdx()
            length=float(np.linalg.norm(positions[mapping[i]]-positions[mapping[j]]))
            radius=sum(table.GetRcovalent(graph.GetAtomWithIdx(a).GetAtomicNum()) for a in (i,j))
            ratios.append(length/radius)
        if not .75<=min(ratios)<=max(ratios)<=1.30:raise ValueError('Implausible covalent bond length in initial XYZ')
        for site,entry in enumerate(topology['para_graft_sites'],1):
            linker=entry['benzyl_atom'];state=case['site_states'][site-1]
            if linker is None:
                if state!='styrene':raise ValueError('Missing declared benzyl handle')
                continue
            atom=graph.GetAtomWithIdx(linker)
            outer=[a for a in atom.GetNeighbors() if a.GetAtomicNum()!=1 and a.GetIdx() not in entry['phenyl_atoms']]
            if len(outer)!=1 or atom.GetTotalNumHs(includeNeighbors=True)!=2:raise ValueError('Benzyl linker must be CH2 attached to exactly one outer atom')
            target=outer[0]
            if state=='grafted':
                if target.GetAtomicNum()!=7 or target.GetFormalCharge()!=1 or not target.GetIsAromatic():raise ValueError('Expected benzyl-to-pyridinium N attachment')
            elif state=='benzyl_chloride':
                if target.GetAtomicNum()!=17 or target.GetDegree()!=1:raise ValueError('An unoccupied handle must retain covalent Cl')
            else:raise ValueError('Unexpected handle state')
            attachments.append(dict(site=site,state=state,linker_atom=linker,outer_atom=target.GetIdx(),bond=graph.GetBondBetweenAtoms(linker,target.GetIdx()).GetIdx()))
        grafts=len(case['occupied_sites']);q=Chem.GetFormalCharge(graph)
        if q!=case['charge'] or q!=2*grafts:raise ValueError('Viologen graph/occupancy charge mismatch')
        salt_mass=Descriptors.MolWt(graph)+grafts*sum(ion['per_graft']*Descriptors.MolWt(Chem.MolFromSmiles(ion['smiles'])) for ion in config['ionic_model']['composition_reference']['counterions'])
        for z,bound in case['composition_bounds'].items():
            capacity=grafts*int(z)*(1.602176634e-19*6.02214076e23)/(3.6*salt_mass)
            if not math.isclose(salt_mass,bound['mass_g_mol'],abs_tol=1e-8) or not math.isclose(capacity,bound['capacity_upper_bound_mAh_g'],abs_tol=1e-8):
                raise ValueError('Graph-derived dry-salt capacity differs from analytic bookkeeping')
        records.append(dict(key=key,charge=q,occupied_sites=case['occupied_sites'],backbone_carbons=topology['backbone_carbon_count'],backbone_edges=topology['backbone_CC_bonds'],chain_ends=topology['chain_ends'],attachments=attachments,minimum_bond_covalent_radius_ratio=min(ratios),maximum_bond_covalent_radius_ratio=max(ratios),selected_MMFF_status=audit['selected']['MMFF_status'],salt_mass_g_mol=salt_mass,all_axis_decision='pending'))
    for cycle in data['cycles']:
        elements=Counter();charge=0
        for key,coefficient in cycle['coefficients'].items():
            if type(coefficient) is not int:raise ValueError('Integer cycle coefficients required')
            for atom in graphs[key].GetAtoms():elements[atom.GetSymbol()]+=coefficient
            charge+=coefficient*Chem.GetFormalCharge(graphs[key])
        if charge or any(elements.values()):raise ValueError('Unbalanced conditional occupancy cycle')
    return dict(status='initial_ionic_graph_XYZ_and_salt_mass_replay_passed',source=str(path.resolve()),source_sha256=sha256_file(path),graphs=len(records),balanced_cycles=len(data['cycles']),records=records,quantum_jobs_submitted=0,native_G_qualified=False,ensemble_sampling_qualified=False,charged_solver_qualified=False,material_admission=False,producer_sha256=sha256_file(__file__))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--preparation',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    from redox.workflows.resources import allocation_guard,host_capacity,cpu_allocation
    with allocation_guard():host_capacity(1,minimum_ram_gib=4);os.sched_setaffinity(0,cpu_allocation(1))
    if a.out.exists():raise FileExistsError('Preserve earlier audits')
    result=check(a.preparation.resolve());atomic_write_json(a.out,result)
    print(json.dumps({k:result[k] for k in ('status','graphs','balanced_cycles','native_G_qualified')}))
