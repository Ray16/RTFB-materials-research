#!/usr/bin/env python3
"""Draw the source-identified precursor and the proposed O-linked local graft."""
import argparse
import json
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import Draw, rdDepictor


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--candidate',type=Path,required=True);p.add_argument('--svg',type=Path,required=True)
    p.add_argument('--png',type=Path);a=p.parse_args()
    row=json.loads(a.candidate.read_text())[0]
    molecules=[Chem.MolFromSmiles(row['route_binding']['precursor_smiles']),Chem.MolFromSmiles(row['product_smiles'])]
    if any(m is None for m in molecules):raise ValueError('Invalid stored molecular graph')
    if Chem.MolToInchiKey(molecules[0])!=row['precursor_inchikey'] or Chem.MolToInchiKey(molecules[1])!=row['product_inchikey']:
        raise ValueError('Graph identity differs')
    highlights=[]
    for i,mol in enumerate(molecules):
        rdDepictor.Compute2DCoords(mol)
        nh=mol.GetSubstructMatches(Chem.MolFromSmarts('[nH]'))
        if len(nh)!=1:raise ValueError('Expected one imide N-H')
        atoms=[nh[0][0]]
        if i==0:
            oh=mol.GetSubstructMatches(Chem.MolFromSmarts('[O;H1]'))
            if len(oh)!=1:raise ValueError('Expected one precursor alcohol')
            atoms.append(oh[0][0])
        highlights.append(atoms)
    labels=['Source-reported precursor: OH and N-H','Proposed O-graft: N-H remains']
    svg=Draw.MolsToGridImage(molecules,molsPerRow=2,subImgSize=(520,330),legends=labels,
                             highlightAtomLists=highlights,useSVG=True)
    a.svg.parent.mkdir(parents=True,exist_ok=True)
    with a.svg.open('x') as stream:stream.write(str(svg))
    if a.png:
        if a.png.exists():raise FileExistsError('Preserve existing figure')
        a.png.parent.mkdir(parents=True,exist_ok=True)
        Draw.MolsToGridImage(molecules,molsPerRow=2,subImgSize=(520,330),legends=labels,
                             highlightAtomLists=highlights).save(a.png)


if __name__=='__main__':main()
