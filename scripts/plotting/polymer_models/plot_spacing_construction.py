#!/usr/bin/env python
"""Draw the exact five-repeat phthalimide occupancy graphs behind the G cycles."""
import argparse,json,os,subprocess,sys
from pathlib import Path
from xml.etree import ElementTree as ET
from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
from redox.core.common import ROOT
from redox.workflows.chunked import file_hash
from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--review',type=Path,default=ROOT/'results/reviews/2026-10-22/spacing_native_complete/validation.json')
    parser.add_argument('--out',type=Path,default=ROOT/'results/figures/polymer_models')
    parser.add_argument('--renderer-python',default=sys.executable,help='Python with CairoSVG; default current interpreter')
    args=parser.parse_args()
    with allocation_guard():host_capacity(1);os.sched_setaffinity(0,cpu_allocation(1))
    data=json.loads(args.review.read_text());rows={tuple(r['occupied_sites']):r for r in data['cases']}
    svg=ET.Element('svg',xmlns='http://www.w3.org/2000/svg',width='1200',height='470',viewBox='0 0 1200 470')
    ET.SubElement(svg,'rect',width='1200',height='470',fill='white')
    proofs=[]
    for index,(occupied,label) in enumerate([((),'a   No grafts'),((2,3),'b   Adjacent grafts (2,3)'),((2,4),'c   One spacer (2,4)')]):
        row=rows[occupied];mol=Chem.MolFromSmiles(row['mapped_smiles']);plain=Chem.Mol(mol)
        backbone={a.GetAtomMapNum()-1000:a.GetIdx() for a in mol.GetAtoms() if 1000<a.GetAtomMapNum()<=1010}
        if len(backbone)!=10 or any(not mol.GetBondBetweenAtoms(backbone[i],backbone[i+1]) for i in range(1,10)):
            raise ValueError('Figure backbone is not the exact open ten-carbon path')
        for a in plain.GetAtoms():a.SetAtomMapNum(0)
        if Chem.MolToSmiles(plain)!=Chem.MolToSmiles(Chem.MolFromSmiles(row['smiles'])):
            raise ValueError('Mapped drawing differs from the computed graph')
        for site in range(1,6):plain.GetAtomWithIdx(backbone[2*site]).SetProp('atomNote',str(site))
        rdDepictor.Compute2DCoords(plain)
        drawer=rdMolDraw2D.MolDraw2DSVG(400,360)
        opts=drawer.drawOptions();opts.bondLineWidth=1.8;opts.annotationFontScale=.75;opts.padding=.08
        atoms=list(backbone.values());bonds=[plain.GetBondBetweenAtoms(backbone[i],backbone[i+1]).GetIdx() for i in range(1,10)]
        drawer.DrawMolecule(plain,highlightAtoms=atoms,highlightBonds=bonds,
            highlightAtomColors={i:(.55,.78,.91) for i in atoms},highlightBondColors={i:(.55,.78,.91) for i in bonds})
        drawer.FinishDrawing();child=ET.fromstring(drawer.GetDrawingText());child.set('x',str(index*400));child.set('y','42')
        svg.append(child)
        text=ET.SubElement(svg,'text',x=str(15+index*400),y='28',style='font-family:sans-serif;font-size:20px;fill:black');text.text=label
        proofs.append(dict(occupied_sites=occupied,smiles=row['smiles'],mapped_smiles=row['mapped_smiles'],site_states=row['site_states'],backbone_atoms=10,backbone_bonds=9,chain_ends=2))
    t=ET.SubElement(svg,'text',x='25',y='445',style='font-family:sans-serif;font-size:17px;fill:black')
    t.text='Sites 1 and 5: plain styrene padding; unoccupied sites 2–4 retain benzyl chloride. Blue: open backbone.'
    args.out.mkdir(parents=True,exist_ok=True);stem=args.out/'spacing_construction';stem.with_suffix('.svg').write_text(ET.tostring(svg,encoding='unicode'))
    subprocess.run([args.renderer_python,'-c',
        'import sys,cairosvg; p=sys.argv[1]; cairosvg.svg2pdf(url=p+".svg",write_to=p+".pdf"); cairosvg.svg2png(url=p+".svg",write_to=p+".png")',str(stem)],check=True)
    stem.with_suffix('.json').write_text(json.dumps(dict(status='exact_computed_graphs_drawn',cases=proofs,
        source_sha256={str(p):file_hash(p) for p in [Path(__file__),args.review]},
        output_sha256={str(stem.with_suffix('.'+suffix)):file_hash(stem.with_suffix('.'+suffix)) for suffix in ('svg','pdf','png')},scope='Molecular graphs, not a sampled 3D conformation'),indent=2)+'\n')


if __name__=='__main__':main()
