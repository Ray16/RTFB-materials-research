"""Deterministic vector structures with explicit model/monomer scope and attachment marks."""
import html
import json
from pathlib import Path

from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D

from candidate_highlights import lead_name


def molecule_svg(smiles,vinyl=False):
    mol=Chem.MolFromSmiles(smiles)
    if mol is None:raise ValueError('Cannot depict an invalid molecular graph')
    query=Chem.MolFromSmarts('[CH2:1]=[CH:2]-c1ccc([CH2:3]-[!#0:4])cc1' if vinyl
                            else '[CH3]c1ccc([CH2:3]-[!#0:4])cc1')
    matches=mol.GetSubstructMatches(query)
    if len(matches)!=1:raise ValueError('Attachment handle is not uniquely identified')
    mapping={atom.GetAtomMapNum():matches[0][atom.GetIdx()] for atom in query.GetAtoms() if atom.GetAtomMapNum()}
    bond=mol.GetBondBetweenAtoms(mapping[3],mapping[4]).GetIdx()
    atoms=[mapping[1],mapping[2]] if vinyl else []
    bonds=[bond];bond_colors={bond:(.95,.65,.30)}
    if vinyl:
        vinyl_bond=mol.GetBondBetweenAtoms(mapping[1],mapping[2]).GetIdx()
        bonds.append(vinyl_bond);bond_colors[vinyl_bond]=(.65,.83,.95)
    drawer=rdMolDraw2D.MolDraw2DSVG(360,250)
    options=drawer.drawOptions();options.fixedFontSize=18;options.bondLineWidth=1.6
    options.padding=.1;options.useBWAtomPalette()
    drawer.DrawMolecule(mol,highlightAtoms=atoms,highlightBonds=bonds,
        highlightAtomColors={a:(.65,.83,.95) for a in atoms},
        highlightBondColors=bond_colors)
    drawer.FinishDrawing()
    return drawer.GetDrawingText().split('?>')[-1]


def structure_plate(rows,report,names,proposal,out,controls=()):
    import cairosvg
    by_id={r['id']:r for r in list(rows)+list(controls)}
    monomers={r['source_capped_id']:r for r in proposal['records']}
    cells=[]
    for i,gid in enumerate(report['highlighted_ids']):
        cells.append(dict(id=gid,smiles=by_id[gid]['exact_grafted_smiles'],
            label=f'{i+1}  '+lead_name(gid,names),scope='computed_capped_local_model',row=0,column=i))
    for i,(gid,label) in enumerate([('design_ndi_ethyl','NDI graft-site model'),('design_pmdi_ethyl','PMDI graft-site model'),('ethylviologen','Viologen graft-site model')]):
        row=by_id[gid]
        cells.append(dict(id=gid,smiles=row['exact_grafted_smiles'],label=f'C{i+1}  '+label,scope='computed_capped_local_model',row=1,column=i))
    parts=['<svg xmlns="http://www.w3.org/2000/svg" width="7.2in" height="5.3in" viewBox="0 0 1080 795">',
           '<rect width="1080" height="795" fill="white"/>',
           '<g font-family="DejaVu Sans,Arial,sans-serif" fill="#202a35">',
           '<text x="20" y="35" font-size="22" font-weight="bold">a  Descriptor leaders · computed capped models</text>',
           '<text x="20" y="420" font-size="22" font-weight="bold">b  Representative graft-site controls · computed capped models</text>']
    for cell in cells:
        x=360*cell['column'];y=65 if cell['row']==0 else 450
        svg=molecule_svg(cell['smiles'],vinyl=False)
        parts.extend([f'<g transform="translate({x},{y})">{svg}</g>',
            f'<text x="{x+180}" y="{y+275}" text-anchor="middle" font-size="18">{html.escape(cell["label"])}</text>'])
    parts.extend(['<text x="20" y="770" font-size="17">Orange bond: attachment to the redox center. These are local models of grafted sites.</text>',
                  '</g></svg>'])
    svg='\n'.join(parts)
    (out/'candidate_structures.svg').write_text(svg)
    cairosvg.svg2pdf(bytestring=svg.encode(),write_to=str(out/'candidate_structures.pdf'))
    cairosvg.svg2png(bytestring=svg.encode(),write_to=str(out/'candidate_structures.png'),dpi=600)
    from PIL import Image
    with Image.open(out/'candidate_structures.png') as raster:
        raster.save(out/'candidate_structures.png',dpi=(600,600))
    (out/'candidate_structures.json').write_text(json.dumps(dict(status='graph_exact_vector_depictions',
        current_route='post-functionalization of a preformed polymer',structures=cells,
        computed_values_transferred_to_vinyl_proposals=False),indent=2)+'\n')
    return monomers
