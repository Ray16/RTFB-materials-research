#!/usr/bin/env python
"""Draw the source-linked TEMPO precursor and capped local product graph."""
import argparse
import json
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import Draw


def draw(route_path,out):
    route=json.loads(route_path.read_text())
    binding=route['evidence']['route_binding']
    precursor=Chem.MolFromSmiles(binding['mapped_precursor_smiles'])
    product=Chem.MolFromSmiles(binding['mapped_product_smiles'])
    maps=binding['new_bond_atom_maps']
    marked=[]
    for molecule,wanted in ((precursor,[maps[1]]),(product,maps)):
        indices=[a.GetIdx() for a in molecule.GetAtoms() if a.GetAtomMapNum() in wanted]
        if len(indices)!=len(wanted):raise ValueError('Mapped attachment atoms missing')
        marked.append(indices)
        for atom in molecule.GetAtoms():atom.SetAtomMapNum(0)
    out.mkdir(parents=True,exist_ok=True)
    options=dict(molsPerRow=2,subImgSize=(560,330),
        legends=['4-Hydroxy-TEMPO precursor','4-Methylbenzyl local graft model'],highlightAtomLists=marked)
    svg=Draw.MolsToGridImage([precursor,product],useSVG=True,**options)
    (out/'tempo_exact_route.svg').write_text(svg)
    image=Draw.MolsToGridImage([precursor,product],useSVG=False,**options)
    image.save(out/'tempo_exact_route.png')
    image.convert('RGB').save(out/'tempo_exact_route.pdf','PDF',resolution=200.0)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--route',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();draw(a.route,a.out)
