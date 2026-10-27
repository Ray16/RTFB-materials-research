#!/usr/bin/env python3
"""Draw a mapped precursor and its proposed local graft with the new bond marked."""
import argparse
import json
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import Draw, rdDepictor


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--candidate',type=Path,required=True)
    p.add_argument('--svg',type=Path,required=True);p.add_argument('--png',type=Path);a=p.parse_args()
    row=json.loads(a.candidate.read_text())[0];binding=row['route_binding']
    molecules=[Chem.MolFromSmiles(binding['mapped_precursor_smiles']),Chem.MolFromSmiles(binding['mapped_product_smiles'])]
    if any(m is None for m in molecules):raise ValueError('Invalid mapped precursor or product')
    product_maps={atom.GetAtomMapNum():atom.GetIdx() for atom in molecules[1].GetAtoms()}
    left,right=binding['new_bond_atom_maps']
    if left not in product_maps or right not in product_maps:raise ValueError('Missing mapped attachment atoms')
    bond=molecules[1].GetBondBetweenAtoms(product_maps[left],product_maps[right])
    if bond is None:raise ValueError('Mapped attachment bond is absent')
    anchors=[[],[product_maps[left],product_maps[right]]];bonds=[[],[bond.GetIdx()]]
    for mol in molecules:
        for atom in mol.GetAtoms():atom.SetAtomMapNum(0)
        rdDepictor.Compute2DCoords(mol)
    if Chem.MolToInchiKey(molecules[0])!=row['precursor_inchikey'] or Chem.MolToInchiKey(molecules[1])!=row['product_inchikey']:
        raise ValueError('Mapped graph does not match declared identity')
    options=dict(molsPerRow=2,subImgSize=(520,330),legends=['Source-identified precursor','Proposed local graft; new bond highlighted'],
                 highlightAtomLists=anchors,highlightBondLists=bonds)
    svg=Draw.MolsToGridImage(molecules,useSVG=True,**options)
    a.svg.parent.mkdir(parents=True,exist_ok=True)
    with a.svg.open('x') as stream:stream.write(str(svg))
    if a.png:
        if a.png.exists():raise FileExistsError('Preserve existing figure')
        a.png.parent.mkdir(parents=True,exist_ok=True);Draw.MolsToGridImage(molecules,**options).save(a.png)


if __name__=='__main__':main()
