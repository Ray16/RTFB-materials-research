#!/usr/bin/env python
"""Localize structural redox relaxation to core and substituent coordinates.

Whole-molecule RMSD alone cannot distinguish a legitimate charge-dependent
conformer preference from a flexible graft's peripheral rotation. This report
keeps all original flags, provides core-aligned displacement and mapped bond
length changes, and never asserts that a polymer freezes a torsion without
an explicit polymer model. No energies are recomputed or ranking flags cleared.
"""
from __future__ import annotations
import argparse
import csv
import json
from pathlib import Path

import numpy as np
from ase.io import read
from rdkit import Chem
from rdkit.Chem import rdDetermineBonds

from redox.core.common import (
    DFT, ROOT, _skeleton, read_manifest, redox_core_xyz_indices,
    smiles_to_xyz_map,
)
from redox.properties.reorg import _couples, lambda_for_couple

OUT=ROOT/'results/validation/physical_reorg_2026-10-10'


def localized_core_indices(smiles,geometry):
    """For carbonyl motifs, select carbonyl-bearing fused rings only.

    The generic helper selects every nontether ring and can therefore count
    an isolated phenyl side chain as redox core. Fused quinone/imide rings and
    their carbonyl oxygens are the specific diagnostic region used here.
    Other families retain the generic ring-based definition.
    """
    generic,tether=redox_core_xyz_indices(smiles,geometry)
    mapping=smiles_to_xyz_map(smiles,geometry)
    if mapping is None:
        return None,tether,'mapping_failed'
    mol=Chem.MolFromSmiles(smiles)
    rings=[set(r) for r in mol.GetRingInfo().AtomRings()]
    carbonyls=[]
    for bond in mol.GetBonds():
        i,j=bond.GetBeginAtomIdx(),bond.GetEndAtomIdx()
        if bond.GetBondType()==Chem.BondType.DOUBLE:
            syms=[mol.GetAtomWithIdx(i).GetSymbol(),mol.GetAtomWithIdx(j).GetSymbol()]
            if sorted(syms)==['C','O']:
                carbonyls.append((i,j) if syms[0]=='C' else (j,i))
    anchors={c for c,o in carbonyls if any(c in r for r in rings)}
    if not anchors:
        return generic,tether,'generic_nontether_rings'
    selected=set(anchors)
    while True:
        extended=selected|set().union(*(r for r in rings if r&selected))
        if extended==selected:
            break
        selected=extended
    selected|={o for c,o in carbonyls if c in selected}
    return sorted(mapping[i] for i in selected),tether,'carbonyl_bearing_fused_rings'


def aligned_coordinates(a,b,idx):
    A=a[idx]; B=b[idx]
    ca,cb=A.mean(0),B.mean(0)
    u,_,vt=np.linalg.svd((B-cb).T@(A-ca))
    d=np.sign(np.linalg.det(u@vt))
    rot=u@np.diag([1,1,d])@vt
    return (b-cb)@rot+ca


def analyze(gid,ox,red,smiles):
    g0=DFT/gid/ox/'opt.xyz'; g1=DFT/gid/red/'opt.xyz'
    aa,bb=read(str(g0)),read(str(g1))
    a=aa.positions[np.array(aa.get_chemical_symbols())!='H']
    b=bb.positions[np.array(bb.get_chemical_symbols())!='H']
    if a.shape!=b.shape:
        raise ValueError('composition mismatch')
    core,tether,core_definition=localized_core_indices(smiles,g0)
    if core is None or len(core)<3:
        raise ValueError('no mapped redox core')
    geometry_mols=[]
    for geom in (g0,g1):
        mol=Chem.MolFromXYZFile(str(geom))
        rdDetermineBonds.DetermineConnectivity(mol)
        geometry_mols.append(_skeleton(Chem.RemoveHs(mol,sanitize=False)))
    matches=geometry_mols[1].GetSubstructMatches(
        geometry_mols[0],uniquify=False,maxMatches=10000)
    if not matches:
        raise ValueError('connectivity mismatch')
    best=None
    for match in matches:
        trial=aligned_coordinates(a,b[list(match)],core)
        delta=np.linalg.norm(trial-a,axis=1)
        core_rmsd=float(np.sqrt(np.mean(delta[core]**2)))
        all_rmsd=float(np.sqrt(np.mean(delta**2)))
        if best is None or (core_rmsd,all_rmsd)<best[:2]:
            best=(core_rmsd,all_rmsd,delta,match)
    cr,ar,delta,match=best
    peripheral=sorted(set(range(len(a)))-set(core))
    rec=dict(id=gid,couple=f'{ox}->{red}',core_rmsd_A=cr,
        core_aligned_all_rmsd_A=ar,
        periphery_rmsd_A=float(np.sqrt(np.mean(delta[peripheral]**2))) if peripheral else None,
        tether_rmsd_A=float(np.sqrt(np.mean(delta[tether]**2))) if tether else None,
        symmetry_maps_checked=len(matches),core_heavy_count=len(core),tether_heavy_count=len(tether),
        core_definition=core_definition,
        interpretation=('small_core_large_periphery' if cr<.15 and ar>.4 else
                        'core_distortion_present' if cr>=.15 else 'small_geometric_change'),
        action='retain_original_flags; conformer energetics and polymer constraints require explicit tests')
    mol=Chem.MolFromSmiles(smiles)
    mapping=smiles_to_xyz_map(smiles,g0)
    lengths=[]
    for bond in mol.GetBonds():
        if bond.GetBondType()!=Chem.BondType.DOUBLE:
            continue
        i,j=bond.GetBeginAtomIdx(),bond.GetEndAtomIdx()
        if sorted([mol.GetAtomWithIdx(i).GetSymbol(),mol.GetAtomWithIdx(j).GetSymbol()])!=['C','O']:
            continue
        ia,ja=mapping[i],mapping[j]
        ro=float(np.linalg.norm(a[ia]-a[ja]))
        rr=float(np.linalg.norm(b[match[ia]]-b[match[ja]]))
        lengths.append(dict(atoms=[ia,ja],oxidized_A=ro,reduced_A=rr,change_A=rr-ro))
    rec['carbonyl_bond_changes']=lengths
    return rec


def main():
    global OUT
    parser=argparse.ArgumentParser()
    parser.add_argument('--out',type=Path,default=OUT)
    args=parser.parse_args();OUT=args.out
    OUT.mkdir(parents=True,exist_ok=True)
    manifest={r['id']:r for r in read_manifest()}
    rows=[]
    for gid,entry in manifest.items():
        try:
            couples=_couples(gid)
        except Exception:
            continue
        for o,r in couples:
            try:
                rec=analyze(gid,o[0],r[0],entry['smiles'])
                raw=lambda_for_couple(gid,o,r)
                rec['original_flag']=raw['flag'] if raw else None
                rows.append(rec)
            except Exception as exc:
                rows.append(dict(id=gid,couple=f'{o[0]}->{r[0]}',error=repr(exc)))
    (OUT/'conformer_localization.json').write_text(json.dumps(rows,indent=2))
    cols=['id','couple','original_flag','core_rmsd_A','core_aligned_all_rmsd_A',
          'periphery_rmsd_A','tether_rmsd_A','symmetry_maps_checked','interpretation','error']
    with (OUT/'conformer_localization.csv').open('w') as f:
        writer=csv.DictWriter(f,fieldnames=cols,extrasaction='ignore')
        writer.writeheader();writer.writerows(rows)
    print(json.dumps(dict(couples=len(rows),errors=sum('error' in r for r in rows),
        flagged_small_core_periphery=[dict(id=r['id'],couple=r['couple'],core=r['core_rmsd_A'],
                periphery=r['periphery_rmsd_A']) for r in rows if r.get('original_flag')
                and r.get('interpretation')=='small_core_large_periphery']),indent=2))


if __name__=='__main__':
    main()
