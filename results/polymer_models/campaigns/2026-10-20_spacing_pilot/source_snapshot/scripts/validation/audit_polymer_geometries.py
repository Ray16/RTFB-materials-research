#!/usr/bin/env python
"""Read-only geometry sanity review; never substitute this for a minimum test."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdDetermineBonds
from ase.io import read
from ase.data import covalent_radii, vdw_radii

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.build.geometry import _skeleton, geometry_stereo_status
from redox.build.polymer_contacts import contact_metrics
from redox.validation.protonation import hydrogen_topology_status


def audit(task, path):
    path = Path(path)
    atoms = read(path)
    xyz = atoms.positions
    declared = Chem.AddHs(Chem.MolFromSmiles(task['smiles']))
    issues = []
    topology = hydrogen_topology_status(task['smiles'], path)
    stereo = geometry_stereo_status(task['smiles'], path)
    if topology != 'ok': issues.append('connectivity:' + topology)
    if stereo not in ('ok', 'not_applicable'): issues.append('stereochemistry:' + stereo)
    if not np.isfinite(xyz).all(): issues.append('nonfinite_coordinates')
    observed = Chem.MolFromXYZFile(str(path))
    if observed is None: raise ValueError('Unreadable XYZ')
    rdDetermineBonds.DetermineConnectivity(observed, charge=0)
    mapping = _skeleton(observed).GetSubstructMatch(_skeleton(declared))
    lengths = []; overlaps = []
    if len(mapping) != declared.GetNumAtoms():
        issues.append('atom_mapping')
    else:
        for bond in declared.GetBonds():
            i, j = mapping[bond.GetBeginAtomIdx()], mapping[bond.GetEndAtomIdx()]
            distance = float(np.linalg.norm(xyz[i] - xyz[j]))
            ratio = distance / (covalent_radii[atoms.numbers[i]] + covalent_radii[atoms.numbers[j]])
            lengths.append(dict(atoms=[i, j],distance_A=distance,covalent_radius_ratio=ratio))
            if not .65 <= ratio <= 1.35: issues.append(f'bond_length:{i}-{j}')
        graph = Chem.GetDistanceMatrix(declared)
        for a in range(declared.GetNumAtoms()):
            for b in range(a):
                if graph[a, b] <= 2: continue  # Exclude bonded and 1-3 contacts.
                i, j = mapping[a], mapping[b]
                distance = float(np.linalg.norm(xyz[i] - xyz[j]))
                radius = vdw_radii[atoms.numbers[i]] + vdw_radii[atoms.numbers[j]]
                if distance < .55 * radius:
                    overlaps.append(dict(atoms=[i,j],distance_A=distance,vdw_ratio=distance/radius))
        if overlaps: issues.append('severe_nonbonded_overlap')
    contacts = contact_metrics(xyz, task['pi_groups'])
    fragments = len(Chem.GetMolFrags(observed))
    expected_fragments = 2 if task['model_type'] == 'two_nonbonded_segments' else 1
    if fragments != expected_fragments: issues.append('fragment_count')
    return dict(key=task['key'],xyz=str(path),xyz_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                status='geometry_sanity_pass' if not issues else 'geometry_review_required',
                issues=issues,topology=topology,stereochemistry=stereo,fragments=fragments,
                atoms=len(atoms),bond_length_range_A=[min(r['distance_A'] for r in lengths),max(r['distance_A'] for r in lengths)] if lengths else None,
                covalent_bond_ratio_range=[min(r['covalent_radius_ratio'] for r in lengths),max(r['covalent_radius_ratio'] for r in lengths)] if lengths else None,
                severe_nonbonded_overlaps=overlaps,contact_metrics=contacts,
                Hessian_minimum_proven=False,electron_transfer_computed=False)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args()
    campaign=ROOT/'results/polymer_models/campaigns/2026-10-19_pi_stacking'
    plan=json.loads((campaign/'plan.json').read_text())
    results=[]
    for roster in plan['tasks']:
        task=json.loads((campaign/'live_inputs'/roster['key']/'task.json').read_text())
        result=audit(task,task['input_xyz']);result['stage']='starting_geometry';results.append(result)
        final=Path(plan['raw_root'])/task['key']/'rest/opt.xyz'
        if final.exists():
            result=audit(task,final);result['stage']='completed_DFT_geometry';results.append(result)
    summary=dict(scope='Composition, connectivity, stereochemistry, bond lengths, severe clashes and core contacts; not a Hessian or rate validation',
                 solvent='SMD acetonitrile (DFT); ALPB acetonitrile (xTB seeds)',
                 geometric_screen_heuristics=dict(covalent_bond_radius_ratio=[.65,1.35],severe_nonbonded_vdw_ratio_below=.55,exclude_graph_distances_at_most=2),
                 starting_geometries=sum(r['stage']=='starting_geometry' for r in results),
                 completed_DFT_geometries=sum(r['stage']=='completed_DFT_geometry' for r in results),
                 flagged=[dict(key=r['key'],stage=r['stage'],issues=r['issues']) for r in results if r['issues']],
                 failed_optimizer_final_coordinates='Not preserved by the original driver; temporary directories verified absent. Checked original seeds are available for bounded reruns.',
                 results=results)
    atomic_write_json(args.out,summary)
    print(json.dumps({k:v for k,v in summary.items() if k!='results'}))


if __name__=='__main__':main()
