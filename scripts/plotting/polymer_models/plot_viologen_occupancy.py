#!/usr/bin/env python3
"""Draw the actual audited charged-polymer graphs, not an inferred triangle."""
import argparse
import json
from pathlib import Path
import re
import subprocess
import sys

from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D

from redox.build.oligomer_topology import chain_topology
from redox.build.polymer_contacts import pi_core_groups
from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json


def draw(source,audit_path,out,renderer):
    data=json.loads(source.read_text());audit=json.loads(audit_path.read_text())
    if audit['status']!='initial_ionic_graph_XYZ_and_salt_mass_replay_passed' or audit['source_sha256']!=sha256_file(source):raise ValueError('The exact preparation must pass its initial-structure audit')
    colors=dict(backbone=(.30,.63,.68),core=(.95,.69,.35),chloride=(.68,.58,.79))
    panels=[];bindings=[]
    for index,(key,title) in enumerate([('all_same__occ_2_3','Adjacent grafts: sites 2 + 3'),('all_same__occ_2_4','One spacer: sites 2 + 4')]):
        row=next(r for r in data['cases'] if r['case']['key']==key);case=row['case'];mol=Chem.MolFromSmiles(case['smiles'])
        topology=chain_topology(case['smiles'],5,[2,3,4]);groups=pi_core_groups(mol,data['fragment'])
        atom_colors={i:colors['backbone'] for i in topology['backbone_atoms']}
        atom_colors.update({i:colors['core'] for group in groups for i in group})
        atom_colors.update({a.GetIdx():colors['chloride'] for a in mol.GetAtoms() if a.GetAtomicNum()==17})
        for site,entry in enumerate(topology['para_graft_sites'],1):mol.GetAtomWithIdx(entry['backbone_atom']).SetProp('atomNote',str(site))
        bond_colors={bond.GetIdx():atom_colors[bond.GetBeginAtomIdx()] for bond in mol.GetBonds() if bond.GetBeginAtomIdx() in atom_colors and bond.GetEndAtomIdx() in atom_colors and atom_colors[bond.GetBeginAtomIdx()]==atom_colors[bond.GetEndAtomIdx()]}
        rdDepictor.Compute2DCoords(mol,nSample=200,nFlipsPerSample=3,sampleSeed=20261031,permuteDeg4Nodes=True)
        drawer=rdMolDraw2D.MolDraw2DSVG(1280,360);options=drawer.drawOptions()
        options.bondLineWidth=2.;options.minFontSize=15;options.maxFontSize=24;options.annotationFontScale=.8;options.padding=.04
        rdMolDraw2D.PrepareAndDrawMolecule(drawer,mol,highlightAtoms=list(atom_colors),highlightBonds=list(bond_colors),highlightAtomColors=atom_colors,highlightBondColors=bond_colors)
        drawer.FinishDrawing();svg=drawer.GetDrawingText();body=re.sub(r'^.*?<svg[^>]*>','',svg,flags=re.S);body=re.sub(r'</svg>\s*$','',body)
        y=120+index*420
        panels.append(f'<text x="40" y="{y-15}" font-size="25" font-weight="bold">{title}</text><g transform="translate(20,{y})">{body}</g>')
        bindings.append(dict(key=key,smiles=case['smiles'],charge=case['charge'],backbone_atoms=topology['backbone_atoms'],redox_core_atoms=groups))
    svg='''<svg xmlns="http://www.w3.org/2000/svg" xmlns:svg="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" width="1320" height="1050" viewBox="0 0 1320 1050"><rect width="1320" height="1050" fill="white"/><g font-family="DejaVu Sans,Arial,sans-serif" fill="#15343b"><text x="40" y="43" font-size="29" font-weight="bold">Exact viologen on an open five-repeat polymer model</text><text x="40" y="77" font-size="19">Same handle positions (2, 3, 4); only occupied sites change. Two grafts give bare charge +4.</text>'''+''.join(panels)+'''<text x="40" y="940" font-size="18">Teal: open carbon backbone  ·  Orange: viologen cores  ·  Purple: residual benzyl chloride</text><text x="40" y="972" font-size="18">Numbers label repeat sites. Terminal repeats 1 and 5 are unfunctionalized styrene padding.</text><text x="40" y="1004" font-size="18">Connectivity drawings; prepared XYZ structures are audited separately. No interaction energy is shown.</text><text x="40" y="1034" font-size="18">Counterions are omitted from coordinates. PF₆⁻ appears only in a separate dry-salt mass convention.</text></g></svg>'''
    if any(out.with_suffix('.'+ext).exists() for ext in ['svg','pdf','png','json']):raise FileExistsError('Preserve earlier model drawings')
    out.parent.mkdir(parents=True,exist_ok=True);out.with_suffix('.svg').write_text(svg)
    program='import pathlib,pymupdf,sys;p=pathlib.Path(sys.argv[1]);d=pymupdf.open(stream=p.with_suffix(".svg").read_bytes(),filetype="svg");pdf=d.convert_to_pdf();p.with_suffix(".pdf").write_bytes(pdf);r=pymupdf.open(stream=pdf,filetype="pdf");r[0].get_pixmap(matrix=pymupdf.Matrix(1,1),alpha=False).save(str(p.with_suffix(".png")))'
    subprocess.run([renderer,'-c',program,str(out)],check=True)
    atomic_write_json(out.with_suffix('.json'),dict(status='actual_audited_ionic_polymer_connectivity_figure',source_sha256=sha256_file(source),audit_sha256=sha256_file(audit_path),producer_sha256=sha256_file(__file__),graphs=bindings,outputs_sha256={ext:sha256_file(out.with_suffix('.'+ext)) for ext in ['svg','pdf','png']},interaction_strength_computed=False,native_G_qualified=False))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('preparation','audit','out'):p.add_argument('--'+key,type=Path,required=True)
    p.add_argument('--renderer-python',default=sys.executable);a=p.parse_args()
    draw(a.preparation.resolve(),a.audit.resolve(),a.out.resolve(),a.renderer_python)
