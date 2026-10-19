"""Illustrative valence/resonance drawings for specified one- and two-electron ladders.

These are chemical schematics, never inferred charge/spin-density predictions.
All transformations preserve atom identity, connectivity and hydrogen count.
"""
from collections import Counter
import html
import itertools
import json
import hashlib
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors


def frozen_kekule(smiles):
    mol=Chem.MolFromSmiles(smiles)
    if mol is None:raise ValueError('Invalid resting-state graph')
    mol=Chem.Mol(mol)
    for atom in mol.GetAtoms():
        count=atom.GetTotalNumHs()
        atom.SetNumExplicitHs(count);atom.SetNoImplicit(True)
    Chem.Kekulize(mol,clearAromaticFlags=True)
    return mol


def alternating_paths(mol,start,end,first,maximum=16):
    """Simple alternating single/double paths, deterministic and bounded."""
    def visit(path,wanted):
        here=path[-1]
        if here==end:
            if len(path)>1:yield path
            return
        if len(path)>maximum:return
        atom=mol.GetAtomWithIdx(here)
        for bond in atom.GetBonds():
            other=bond.GetOtherAtomIdx(here)
            if other in path or int(bond.GetBondTypeAsDouble())!=wanted:continue
            if mol.GetAtomWithIdx(other).GetAtomicNum()!=6 and other!=end:continue
            yield from visit(path+[other],3-wanted)
    yield from visit([start],first)


def flip_path(mol,path):
    for first,second in zip(path,path[1:]):
        bond=mol.GetBondBetweenAtoms(first,second)
        bond.SetBondType(Chem.BondType.DOUBLE if bond.GetBondType()==Chem.BondType.SINGLE else Chem.BondType.SINGLE)


def validated_state(original,mol,charge,radicals):
    Chem.SanitizeMol(mol)
    if Chem.GetFormalCharge(mol)!=charge:raise ValueError('Incorrect molecular charge')
    if sum(a.GetNumRadicalElectrons() for a in mol.GetAtoms())!=radicals:raise ValueError('Incorrect radical electron count')
    atoms=lambda m:Counter(a.GetAtomicNum() for a in Chem.AddHs(m).GetAtoms())
    if atoms(mol)!=atoms(original):raise ValueError('Redox drawing changed elemental/H composition')
    edges=lambda m:{tuple(sorted((b.GetBeginAtomIdx(),b.GetEndAtomIdx()))) for b in m.GetBonds()}
    if edges(mol)!=edges(original):raise ValueError('Redox drawing changed molecular connectivity')
    return mol


def reduction_forms(smiles):
    """Return resting/radical/reduced illustrative forms; refuse unknown chemistry."""
    original=Chem.MolFromSmiles(smiles);base=frozen_kekule(smiles);q=Chem.GetFormalCharge(base)
    if q==2:
        ends=[a.GetIdx() for a in base.GetAtoms() if a.GetAtomicNum()==7 and a.GetFormalCharge()==1]
        if len(ends)!=2:raise ValueError('Dication is not the supported viologen motif')
        paths=[p for p in alternating_paths(base,ends[0],ends[1],2) if len(p)%2==0]
        for path in sorted(paths,key=lambda p:(len(p),p)):
            try:
                reduced=Chem.Mol(base);flip_path(reduced,path)
                for index in ends:reduced.GetAtomWithIdx(index).SetFormalCharge(0)
                validated_state(original,reduced,0,0)
                # An odd alternating prefix produces a carbon-centered formal
                # resonance radical without asserting its physical localization.
                length=(len(path)-1)//2
                if length%2==0:length-=1
                radical=Chem.Mol(base);flip_path(radical,path[:length+1])
                radical.GetAtomWithIdx(path[0]).SetFormalCharge(0)
                radical.GetAtomWithIdx(path[length]).SetNumRadicalElectrons(1)
                validated_state(original,radical,1,1)
                return [original,radical,reduced]
            except ValueError:continue
        raise ValueError('No valid viologen resonance path')
    if q!=0:raise ValueError('Unsupported resting charge for a formal resonance drawing')
    carbonyls=[]
    for bond in base.GetBonds():
        a,b=bond.GetBeginAtom(),bond.GetEndAtom()
        if bond.GetBondType()!=Chem.BondType.DOUBLE:continue
        if a.GetAtomicNum()==8 and b.GetAtomicNum()==6:carbonyls.append((b.GetIdx(),a.GetIdx()))
        if b.GetAtomicNum()==8 and a.GetAtomicNum()==6:carbonyls.append((a.GetIdx(),b.GetIdx()))
    for (ca,oa),(cb,ob) in itertools.combinations(sorted(carbonyls),2):
        na={n.GetIdx() for n in base.GetAtomWithIdx(ca).GetNeighbors() if n.GetAtomicNum()==7}
        nb={n.GetIdx() for n in base.GetAtomWithIdx(cb).GetNeighbors() if n.GetAtomicNum()==7}
        if na&nb:continue  # one carbonyl from each imide unit when applicable
        paths=[p for p in alternating_paths(base,ca,cb,1) if len(p)%2==0]
        for path in sorted(paths,key=lambda p:(len(p),p)):
            try:
                reduced=Chem.Mol(base)
                for carbon,oxygen in ((ca,oa),(cb,ob)):
                    reduced.GetBondBetweenAtoms(carbon,oxygen).SetBondType(Chem.BondType.SINGLE)
                    reduced.GetAtomWithIdx(oxygen).SetFormalCharge(-1)
                flip_path(reduced,path);validated_state(original,reduced,-2,0)
                radical=Chem.Mol(base)
                radical.GetBondBetweenAtoms(ca,oa).SetBondType(Chem.BondType.SINGLE)
                radical.GetAtomWithIdx(oa).SetFormalCharge(-1)
                radical.GetAtomWithIdx(ca).SetNumRadicalElectrons(1)
                validated_state(original,radical,-1,1)
                return [original,radical,reduced]
            except ValueError:continue
    raise ValueError('No supported valid carbonyl reduction drawing')


def write_schemes(rows,report,names,proposal,out,root):
    import cairosvg
    from PIL import Image
    from pypdf import PdfReader,PdfWriter
    from rdkit.Chem.Draw import rdMolDraw2D
    from candidate_highlights import lead_name
    by_id={r['id']:r for r in rows}
    computed=[dict(by_id[gid],smiles=by_id[gid]['exact_grafted_smiles'],
                   drawing_name=f'{i+1}  '+lead_name(gid,names)) for i,gid in enumerate(report['highlighted_ids'])]
    pilot=sorted([r for r in proposal['records'] if r['initial_modeling_pilot']],key=lambda r:r['source_capped_id'])
    pilot_names={'design_ndi_ethyl':'NDI vinyl monomer','design_pmdi_ethyl':'PMDI vinyl monomer','ethylviologen':'Viologen vinyl monomer'}
    for i,row in enumerate(pilot):row=dict(row,drawing_name=f'P{i+1}  '+pilot_names[row['source_capped_id']]);pilot[i]=row
    records=[];pdfs=[]
    for prefix,group,is_computed in [('redox_states_computed',computed,True),('redox_states_vinyl',pilot,False)]:
        parts=['<svg xmlns="http://www.w3.org/2000/svg" width="7.2in" height="7in" viewBox="0 0 1080 1050">',
            '<defs><marker id="arrow" markerWidth="8" markerHeight="8" refX="7" refY="4" orient="auto"><path d="M0,0 L8,4 L0,8" fill="#596573"/></marker></defs>',
            '<rect width="1080" height="1050" fill="white"/>',
            '<g font-family="DejaVu Sans,Arial,sans-serif" fill="#202a35">',
            '<text x="25" y="30" font-size="22" font-weight="bold">'+('Computed capped-model charge/spin ladders' if is_computed else 'Proposed vinyl-monomer ladders · uncomputed')+'</text>']
        for i,row in enumerate(group):
            forms=reduction_forms(row['smiles']);states=row['states']
            if isinstance(states,str):states=json.loads(states)
            states=sorted(states,key=lambda state:abs(state[3]))
            if len(states)!=3:raise ValueError('State diagram requires the declared two-electron ladder')
            top=65+310*i
            parts.append(f'<text x="25" y="{top}" font-size="20" font-weight="bold">{html.escape(row["drawing_name"])}</text>')
            for j,(mol,state) in enumerate(zip(forms,states)):
                label,q,m,ne=state;evidence={}
                if is_computed:
                    path=root/'calcs/dft'/row['id']/label/'result.json'
                    raw=json.loads(path.read_text());q,m=raw['charge'],raw['mult']
                    evidence[str(path)]=hashlib.sha256(path.read_bytes()).hexdigest()
                validated_state(Chem.MolFromSmiles(row['smiles']),mol,q,m-1)
                drawer=rdMolDraw2D.MolDraw2DSVG(290,215);options=drawer.drawOptions()
                options.fixedFontSize=17;options.bondLineWidth=1.5;options.padding=.09
                drawer.DrawMolecule(mol);drawer.FinishDrawing();svg=drawer.GetDrawingText().split('?>')[-1]
                x=30+365*j;y=top+15
                parts.append(f'<g transform="translate({x},{y})">{svg}</g>')
                state_name=(['Dication','Radical cation','Neutral'] if Chem.GetFormalCharge(forms[0])==2
                            else ['Neutral','Radical anion','Dianion'])[j]
                text=state_name+(f' · q = {q:+d}' if q else ' · q = 0')
                parts.append(f'<text x="{x+145}" y="{top+250}" text-anchor="middle" font-size="17">{html.escape(text)}</text>')
                parts.append(f'<text x="{x+145}" y="{top+273}" text-anchor="middle" font-size="16">M = {m}'+(' (selected model)' if is_computed else ' (proposed)')+'</text>')
                if j<2:
                    ax=x+309;ay=top+122
                    parts.append(f'<line x1="{ax}" y1="{ay}" x2="{ax+35}" y2="{ay}" stroke="#596573" stroke-width="1.5" marker-end="url(#arrow)"/>')
                    parts.append(f'<text x="{ax+17}" y="{ay-15}" text-anchor="middle" font-size="17">+ e−</text>')
                    value=row.get('E1_vs_Fc_V' if j==0 else 'E2_vs_Fc_V')
                    potential=f'{value:.3f} V' if is_computed and value is not None else 'pending'
                    parts.append(f'<text x="{ax+17}" y="{ay+32}" text-anchor="middle" font-size="16">{potential}</text>')
                records.append(dict(id=row['id'],state=label,charge=q,mult=m,n_e=ne,
                    depicted_smiles=Chem.MolToSmiles(mol,isomericSmiles=True),source_smiles=row['smiles'],
                    radicals=sum(a.GetNumRadicalElectrons() for a in mol.GetAtoms()),
                    scope='computed_selected_spin_ladder' if is_computed else 'proposed_uncomputed_ladder',
                    depiction='one illustrative formal resonance representation, not calculated density',
                    evidence_sha256=evidence))
        parts.extend(['<text x="25" y="1017" font-size="16">Illustrative resonance forms; charge/spin may be delocalized. No proton transfer is implied.</text>',
            '<text x="25" y="1040" font-size="16">'+('Potentials: calculated MeCN vs Fc/Fc+. Arrows show electron count, not demonstrated reversibility.' if is_computed else 'Vinyl-state properties and polymer performance are uncomputed. M is the proposed multiplicity.')+'</text>',
            '</g></svg>'])
        svg='\n'.join(parts);(out/f'{prefix}.svg').write_text(svg)
        cairosvg.svg2pdf(bytestring=svg.encode(),write_to=str(out/f'{prefix}.pdf'))
        cairosvg.svg2png(bytestring=svg.encode(),write_to=str(out/f'{prefix}.png'),dpi=600)
        with Image.open(out/f'{prefix}.png') as raster:raster.save(out/f'{prefix}.png',dpi=(600,600))
        pdfs.append(out/f'{prefix}.pdf')
    writer=PdfWriter()
    for path in pdfs:writer.add_page(PdfReader(path).pages[0])
    with (out/'redox_states.pdf').open('wb') as stream:writer.write(stream)
    result=dict(status='validated_formal_redox_state_schematics',states=len(records),records=records,
        computed_spin_density=False,formal_charge_placement_is_not_QM_localization=True,
        same_atoms_H_count_and_connectivity_verified=True,polymerization_validation=False,
        references=['https://pubs.rsc.org/en/content/articlehtml/2024/sc/d4sc02431a',
                    'https://pubs.acs.org/doi/10.1021/ja303173n'])
    (out/'redox_states.json').write_text(json.dumps(result,indent=2)+'\n')
    return result
