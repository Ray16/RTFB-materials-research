#!/usr/bin/env python3
"""Draw the verified attachment graph, with no inferred charge/spin localization."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys

from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D


def draw(source,route_path,output,renderer):
    candidate=json.loads(source.read_text());route=json.loads(route_path.read_text())
    if route['route_status']!='proceed' or not route['exact_route_binding_verified']:
        raise ValueError('An exact mapped attachment pass is required')
    binding=route['evidence']['route_binding'];molecule=Chem.MolFromSmiles(binding['mapped_product_smiles'])
    maps={atom.GetAtomMapNum():atom.GetIdx() for atom in molecule.GetAtoms()}
    atoms=[maps[n] for n in binding['new_bond_atom_maps']]
    bond=molecule.GetBondBetweenAtoms(*atoms)
    if bond is None:raise ValueError('Mapped attachment bond missing')
    for atom in molecule.GetAtoms():atom.SetAtomMapNum(0)
    if (Chem.MolToSmiles(molecule)!=candidate['product_smiles']
            or Chem.MolToInchiKey(molecule)!=candidate['product_inchikey']):
        raise ValueError('Drawing and checked candidate identities differ')
    drawer=rdMolDraw2D.MolDraw2DSVG(900,360);options=drawer.drawOptions()
    options.bondLineWidth=2.;options.minFontSize=20;options.maxFontSize=26
    rdMolDraw2D.PrepareAndDrawMolecule(drawer,molecule,highlightAtoms=atoms,highlightBonds=[bond.GetIdx()],
        highlightAtomColors={i:(.86,.55,.30) for i in atoms},highlightBondColors={bond.GetIdx():(.86,.55,.30)})
    drawer.FinishDrawing();output.parent.mkdir(parents=True,exist_ok=True)
    svg=output.with_suffix('.svg');svg.write_text(drawer.GetDrawingText())
    program='import pathlib,pymupdf,sys;p=pathlib.Path(sys.argv[1]);d=pymupdf.open(stream=p.read_bytes(),filetype="svg");pathlib.Path(sys.argv[2]).write_bytes(d.convert_to_pdf())'
    subprocess.run([renderer,'-c',program,str(svg),str(output)],check=True)
    report=dict(source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),route_sha256=hashlib.sha256(route_path.read_bytes()).hexdigest(),
        product_inchikey=candidate['product_inchikey'],new_bond_atom_maps=binding['new_bond_atom_maps'],
        figure_sha256=hashlib.sha256(output.read_bytes()).hexdigest(),renderer_python=renderer,
        scope='Two-dimensional resting-state attachment connectivity; no orbital/spin localization or native 3D geometry inferred')
    output.with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n');return report


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('source','route','out'):p.add_argument('--'+key,type=Path,required=True)
    p.add_argument('--renderer-python',default=sys.executable);a=p.parse_args()
    print(json.dumps(draw(a.source,a.route,a.out,a.renderer_python)))
