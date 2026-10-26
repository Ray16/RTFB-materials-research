#!/usr/bin/env python
"""Prepare fixed-parent local motif pairs with explicit benzyl-to-methyl caps.

One declared wildcard denotes the primary benzyl attachment carbon. The complete
motif and that carbon are retained, including their original hydrogens. Only the
bond to the polymer aryl group is replaced by C-H. No core is moved or relaxed.
"""
import argparse
import json
import os
from pathlib import Path
import time

import numpy as np
from rdkit import Chem

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.build.polymer_contacts import pi_core_groups,contact_metrics
from redox.workflows import polymer_campaign as api


def fragments(parent,template,cap_length_A=1.09):
    if not np.isfinite(cap_length_A) or not 1.0<=cap_length_A<=1.2:
        raise ValueError('A declared C-H cap length between 1.0 and 1.2 Angstrom is required')
    query=Chem.MolFromSmarts(template)
    if query is None:raise ValueError('Invalid motif template')
    wild=[a.GetIdx() for a in query.GetAtoms() if a.GetAtomicNum()==0]
    if len(wild)!=1 or query.GetAtomWithIdx(wild[0]).GetDegree()!=1:
        raise ValueError('Exactly one terminal attachment wildcard is supported')
    coords=np.asarray(parent.GetConformer().GetPositions())
    if coords.shape!=(parent.GetNumAtoms(),3) or not np.isfinite(coords).all():raise ValueError('Finite parent coordinates required')
    matches={}
    for match in parent.GetSubstructMatches(query,uniquify=False,maxMatches=10000):
        matches.setdefault(frozenset(match),match)
    groups=sorted(matches.values(),key=lambda m:min(m));result=[];used=set()
    for match in groups:
        anchor=match[wild[0]];atom=parent.GetAtomWithIdx(anchor);heavy=set(match)
        if used&heavy:raise ValueError('Overlapping motif assignments cannot be fragmented independently')
        used|=heavy
        outside=[(i,n.GetIdx()) for i in heavy for n in parent.GetAtomWithIdx(i).GetNeighbors()
                 if n.GetAtomicNum()!=1 and n.GetIdx() not in heavy]
        if (atom.GetAtomicNum()!=6 or atom.GetIsAromatic() or atom.GetHybridization()!=Chem.HybridizationType.SP3
                or sum(n.GetAtomicNum()==1 for n in atom.GetNeighbors())!=2
                or len(outside)!=1 or outside[0][0]!=anchor
                or not parent.GetAtomWithIdx(outside[0][1]).GetIsAromatic()
                or parent.GetBondBetweenAtoms(*outside[0]).GetBondType()!=Chem.BondType.SINGLE):
            raise ValueError('Each motif must have one primary benzyl-to-aryl cut bond and two anchor hydrogens')
        keep=heavy|{n.GetIdx() for i in heavy for n in parent.GetAtomWithIdx(i).GetNeighbors() if n.GetAtomicNum()==1}
        indices=sorted(keep);mapping={old:new for new,old in enumerate(indices)}
        target=Chem.RWMol(parent)
        for i in reversed(range(parent.GetNumAtoms())):
            if i not in keep:target.RemoveAtom(i)
        cap=target.AddAtom(Chem.Atom(1));target.AddBond(mapping[anchor],cap,Chem.BondType.SINGLE)
        target=target.GetMol();Chem.SanitizeMol(target)
        conf=Chem.Conformer(target.GetNumAtoms())
        for old,new in mapping.items():conf.SetAtomPosition(new,coords[old])
        direction=coords[outside[0][1]]-coords[anchor];length=np.linalg.norm(direction)
        if not .9<length<2.:raise ValueError('Parent benzyl-aryl bond length is implausible')
        conf.SetAtomPosition(cap,coords[anchor]+direction*cap_length_A/length)
        target.RemoveAllConformers();target.AddConformer(conf)
        if len(Chem.GetMolFrags(target))!=1 or target.GetNumAtoms()!=len(indices)+1:
            raise ValueError('Capped motif is not one connected component')
        charge=sum(parent.GetAtomWithIdx(i).GetFormalCharge() for i in keep)
        if Chem.GetFormalCharge(target)!=charge:raise ValueError('Capping changed the motif formal charge')
        result.append(dict(mol=target,parent_indices=indices,parent_to_fragment=mapping,
            anchor_parent_index=anchor,cut_parent_bond=list(outside[0]),cap_index=cap,
            cap_length_A=float(cap_length_A),charge=charge,core_parent_indices=sorted(heavy-{anchor})))
    if not result:raise ValueError('No complete grafted motif matches the attachment template')
    return result


def xyz(mol):
    coords=mol.GetConformer().GetPositions()
    return str(mol.GetNumAtoms())+'\nFrozen parent motif; constructed cap only\n'+''.join(
        f'{a.GetSymbol()} {r[0]:.12f} {r[1]:.12f} {r[2]:.12f}\n' for a,r in zip(mol.GetAtoms(),coords))


def prepare(source,out,cap_lengths=(1.07,1.09,1.11)):
    data=json.loads(source.read_text());api.verify_sources(data['source_sha256'])
    if len(data['states'])!=4 or data['conformer_convergence_established']:
        raise ValueError('The four reviewed finite-control frames are required')
    out.mkdir(parents=True,exist_ok=False);entries=[];sources={str(p):api.digest(p) for p in
        [source,Path(__file__).resolve(),ROOT/'src/redox/build/polymer_contacts.py']}
    for row in data['states']:
        parent=Chem.AddHs(Chem.MolFromSmiles(row['task']['smiles']))
        if [a.GetAtomicNum() for a in parent.GetAtoms()]!=row['atomic_numbers']:
            raise ValueError('Parent graph and cached atom order differ')
        coords=np.asarray(row['coordinates_A']);conf=Chem.Conformer(len(coords))
        for i,r in enumerate(coords):conf.SetAtomPosition(i,r)
        parent.AddConformer(conf);groups=pi_core_groups(parent,row['task']['fragment'])
        if len(groups)!=2 or {frozenset(g) for g in groups}!={frozenset(g) for g in row['groups']}:
            raise ValueError('The two parent pi-core assignments differ')
        measured=contact_metrics(coords,groups)[0]
        for key in ('centroid_distance_A','interplane_angle_deg','plane_separation_A'):
            if abs(measured[key]-row['contact'][key])>1e-8:raise ValueError('Parent contact geometry differs')
        for cap_length in cap_lengths:
            parts=fragments(parent,row['task']['fragment'],cap_length)
            if len(parts)!=2:raise ValueError('Exactly two independent motifs are required')
            pair=Chem.CombineMols(parts[0]['mol'],parts[1]['mol']);molecules=[parts[0]['mol'],parts[1]['mol'],pair]
            if len(Chem.GetMolFrags(pair))!=2:raise ValueError('Local pair acquired an intersite covalent bond')
            offset=parts[0]['mol'].GetNumAtoms();mapped=[]
            for part,shift in zip(parts,(0,offset)):
                source_group=next(g for g in groups if set(g)<=set(part['parent_indices']))
                mapped.append([part['parent_to_fragment'][i]+shift for i in source_group])
            contact=contact_metrics(pair.GetConformer().GetPositions(),mapped)[0]
            if any(abs(contact[k]-measured[k])>1e-8 for k in ('centroid_distance_A','interplane_angle_deg','plane_separation_A')):
                raise ValueError('Capping moved a motif core')
            key=row['key']+f'__cap{round(cap_length*100):03d}';folder=out/key;folder.mkdir()
            rows=[]
            for label,mol in zip(('A','B','AB'),molecules):
                path=folder/(label+'.xyz');path.write_text(xyz(mol))
                sdf=folder/(label+'.mol');sdf.write_text(Chem.MolToMolBlock(mol,kekulize=False))
                rows.append(dict(id=label,atoms=mol.GetNumAtoms(),charge=Chem.GetFormalCharge(mol),
                    smiles=Chem.MolToSmiles(Chem.RemoveHs(mol),isomericSmiles=True),
                    geometry=str(path),geometry_sha256=api.digest(path),molecular_graph=str(sdf),graph_sha256=api.digest(sdf)))
            provenance=[{k:v for k,v in part.items() if k not in ('mol','parent_to_fragment')} for part in parts]
            entries.append(dict(key=key,parent_key=row['key'],parent_task=row['task'],
                cap_length_A=float(cap_length),fragments=provenance,molecules=rows,contact=contact,
                energies_eV=None,interaction_energy_kJ_mol=None,free_energy_computed=False))
    manifest=dict(status='fixed_parent_local_motif_controls_prepared_no_quantum_values',prepared=time.time(),
        source_sha256=sources,entries=entries,frames=4,cap_lengths_A=list(cap_lengths),
        electronic_tasks=3*len(entries),quantum_jobs_submitted_at_preparation=0,
        interaction_definition='E(AB at parent relative pose)-E(A frozen)-E(B frozen), matched electronic method and medium',
        cap_model='Primary benzyl carbon retained as methyl; aryl bond replaced by H along original bond direction',
        source_frame='Original graph-ordered parent coordinates; no rotations, translations or relaxation',
        planned_media=['gas','acetonitrile'],basis_superposition_correction_computed=False,
        primary_candidate_admission=False,electronic_coupling_computed=False,conformer_population_inferred=False,
        scope='Frozen methyl-capped local motif interaction control; differs from full polymer nonadditivity and free energy')
    atomic_write_json(out/'manifest.json',manifest);return manifest


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity
    with allocation_guard():host_capacity(1);os.sched_setaffinity(0,cpu_allocation(1))
    row=prepare(a.source.resolve(),a.out.resolve());print(json.dumps(dict(status=row['status'],frames=row['frames'],electronic_tasks=row['electronic_tasks'])))


if __name__=='__main__':main()
