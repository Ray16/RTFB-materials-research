#!/usr/bin/env python
"""Independently audit local caps, bond cuts and exact frozen XYZ coordinates."""
import argparse
import json
import os
from pathlib import Path
import time

import numpy as np
from rdkit import Chem

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api


def bonds(mol,offset=0):
    return {(min(b.GetBeginAtomIdx(),b.GetEndAtomIdx())+offset,
             max(b.GetBeginAtomIdx(),b.GetEndAtomIdx())+offset,str(b.GetBondType())) for b in mol.GetBonds()}


def review(campaign,out):
    path=campaign/'manifest.json';manifest=json.loads(path.read_text());api.verify_sources(manifest['source_sha256'])
    source=next(Path(p) for p in manifest['source_sha256'] if Path(p).name=='spacing_conformers.json')
    parents={row['key']:row for row in json.loads(source.read_text())['states']};checked=[]
    if (manifest['primary_candidate_admission'] or manifest['electronic_coupling_computed']
            or manifest['quantum_jobs_submitted_at_preparation'] or len(manifest['entries'])!=12):
        raise ValueError('Expected preparation-only four-frame/three-cap controls')
    proof={**manifest['source_sha256'],str(path):api.digest(path),str(Path(__file__).resolve()):api.digest(__file__)}
    for row in manifest['entries']:
        if row['energies_eV'] is not None or row['interaction_energy_kJ_mol'] is not None or row['free_energy_computed']:
            raise ValueError('Prepared structural controls cannot supply energies')
        parent=parents[row['parent_key']];mol=Chem.AddHs(Chem.MolFromSmiles(parent['task']['smiles']))
        xyz=np.asarray(parent['coordinates_A']);molecules=[];positions=[]
        for item in row['molecules']:
            geometry=Path(item['geometry']);graph=Path(item['molecular_graph'])
            if api.digest(geometry)!=item['geometry_sha256'] or api.digest(graph)!=item['graph_sha256']:
                raise ValueError('Prepared geometry or graph bytes changed')
            proof[str(geometry)]=api.digest(geometry);proof[str(graph)]=api.digest(graph)
            lines=geometry.read_text().splitlines();atoms=[line.split() for line in lines[2:]]
            molecule=Chem.MolFromMolFile(str(graph),removeHs=False)
            if molecule is None or len(atoms)!=int(lines[0]) or len(atoms)!=item['atoms']:
                raise ValueError('Graph or XYZ atom coverage differs')
            coords=np.asarray([[float(v) for v in line[1:]] for line in atoms])
            if coords.shape!=(molecule.GetNumAtoms(),3) or not np.isfinite(coords).all():raise ValueError('Invalid XYZ geometry')
            if [a.GetSymbol() for a in molecule.GetAtoms()]!=[a[0] for a in atoms] or Chem.GetFormalCharge(molecule)!=item['charge']:
                raise ValueError('XYZ/graph elements or formal charges differ')
            molecules.append(molecule);positions.append(coords)
        for part,fragment,coords in zip(row['fragments'],molecules[:2],positions[:2]):
            indices=part['parent_indices'];mapping={old:new for new,old in enumerate(indices)}
            anchor,external=part['cut_parent_bond'];cap=part['cap_index']
            if (indices!=sorted(set(indices)) or cap!=len(indices) or len(coords)!=len(indices)+1
                    or not np.allclose(coords[:-1],xyz[indices],atol=1e-10,rtol=0)):
                raise ValueError('A retained motif or hydrogen moved from its parent')
            expected={(min(mapping[b.GetBeginAtomIdx()],mapping[b.GetEndAtomIdx()]),
                       max(mapping[b.GetBeginAtomIdx()],mapping[b.GetEndAtomIdx()]),str(b.GetBondType()))
                for b in mol.GetBonds() if b.GetBeginAtomIdx() in mapping and b.GetEndAtomIdx() in mapping}
            expected.add((mapping[anchor],cap,'SINGLE'))
            if bonds(fragment)!=expected or fragment.GetAtomWithIdx(cap).GetAtomicNum()!=1:
                raise ValueError('A retained bond or explicit cap bond differs')
            direction=xyz[external]-xyz[anchor];target=xyz[anchor]+direction/np.linalg.norm(direction)*row['cap_length_A']
            if not np.allclose(coords[cap],target,atol=1e-10,rtol=0):raise ValueError('Cap length or original bond direction differs')
            if not mol.GetAtomWithIdx(external).GetIsAromatic() or Chem.GetMolFrags(fragment)!=(tuple(range(len(coords))),):
                raise ValueError('Invalid backbone cut or disconnected motif')
        a,b,pair=molecules;x,y,z=positions;offset=a.GetNumAtoms()
        if (not np.array_equal(z,np.concatenate([x,y])) or len(Chem.GetMolFrags(pair))!=2
                or bonds(pair)!=(bonds(a)|bonds(b,offset))):
            raise ValueError('The local pair changed pose or acquired an intersite bond')
        minimum=float(np.min(np.linalg.norm(x[:,None,:]-y[None,:,:],axis=2)))
        if minimum<.75:raise ValueError('Constructed local fragments have a severe interfragment clash')
        checked.append(dict(key=row['key'],fragment_atoms=[a.GetNumAtoms(),b.GetNumAtoms()],pair_atoms=pair.GetNumAtoms(),
            cap_length_A=row['cap_length_A'],minimum_interfragment_distance_A=minimum,parent_coordinates_preserved=True))
    if ({(r['parent_key'],r['cap_length_A']) for r in manifest['entries']}
            !={(key,length) for key in parents for length in (1.07,1.09,1.11)}):
        raise ValueError('Missing parent frame or declared cap sensitivity')
    result=dict(status='fixed_parent_local_motif_construction_independently_checked',checked=time.time(),
        frames=4,cap_lengths_A=[1.07,1.09,1.11],molecular_inputs=36,models=checked,source_sha256=proof,
        geometry_of_record='Full precision XYZ; MOL files provide bond topology and rounded drawing coordinates',
        electronic_energies_computed=False,interaction_strength_established=False,primary_candidate_admission=False,
        electronic_coupling_computed=False,scope='Structural construction and pose audit only; no energetic or population inference')
    out.parent.mkdir(parents=True,exist_ok=True);atomic_write_json(out,result);return result


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--campaign',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity
    with allocation_guard():host_capacity(1);os.sched_setaffinity(0,cpu_allocation(1))
    result=review(a.campaign.resolve(),a.out.resolve());print(json.dumps({k:result[k] for k in ['status','frames','molecular_inputs']}))


if __name__=='__main__':main()
