#!/usr/bin/env python3
"""Draw exact product graphs in the anthraquinone benzyl-ester matrix."""
import argparse
import json
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import Draw,rdDepictor


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--matrix',type=Path,required=True)
    p.add_argument('--svg',type=Path,required=True);p.add_argument('--png',type=Path);a=p.parse_args()
    rows=json.loads(a.matrix.read_text())['candidates']
    if len(rows)!=4:raise ValueError('Expected four exact products')
    mols=[];labels=[]
    for row in rows:
        molecule=Chem.MolFromSmiles(row['product_smiles'])
        if molecule is None or Chem.MolToInchiKey(molecule)!=row['product_inchikey']:raise ValueError('Product identity differs')
        rdDepictor.Compute2DCoords(molecule);mols.append(molecule);f=row['design_factors']
        labels.append(f"{f['linker']} carboxylate; {f['substituent']}")
    options=dict(molsPerRow=2,subImgSize=(500,330),legends=labels)
    svg=Draw.MolsToGridImage(mols,useSVG=True,**options);a.svg.parent.mkdir(parents=True,exist_ok=True)
    with a.svg.open('x') as stream:stream.write(str(svg))
    if a.png:
        if a.png.exists():raise FileExistsError('Preserve existing figure')
        a.png.parent.mkdir(parents=True,exist_ok=True);Draw.MolsToGridImage(mols,**options).save(a.png)
