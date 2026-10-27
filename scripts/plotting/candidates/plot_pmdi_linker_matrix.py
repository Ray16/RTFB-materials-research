#!/usr/bin/env python3
"""Draw the exact local product graphs in the proposed PMDI comparison."""
import argparse
import json
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import Draw, rdDepictor


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--matrix',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--png',type=Path)
    args=p.parse_args();rows=json.loads(args.matrix.read_text())['candidates']
    if len(rows)!=4:raise ValueError('Expected four exact product graphs')
    mols=[];labels=[]
    for row in rows:
        molecule=Chem.MolFromSmiles(row['product_smiles'])
        if molecule is None or Chem.MolToInchiKey(molecule)!=row['product_inchikey']:
            raise ValueError('Product graph differs from matrix identity')
        rdDepictor.Compute2DCoords(molecule);mols.append(molecule)
        factors=row['design_factors'];origin='existing' if row['source_library_id'] else 'new'
        labels.append(f"{factors['linker']} linker; edge {factors['edge_substituent']} ({origin})")
    image=Draw.MolsToGridImage(mols,molsPerRow=2,subImgSize=(500,320),legends=labels,useSVG=True)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    with args.out.open('x') as stream:stream.write(str(image))
    if args.png:
        if args.png.exists():raise FileExistsError('Preserve existing figure')
        args.png.parent.mkdir(parents=True,exist_ok=True)
        Draw.MolsToGridImage(mols,molsPerRow=2,subImgSize=(500,320),legends=labels).save(args.png)


if __name__=='__main__':
    main()
