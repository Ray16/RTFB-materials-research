"""One maintained, graph-exact NDI/cata-NDI comparison for the research story."""
import html,json,math
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from candidate_structure_plate import molecule_svg

IDS=('design_ndi_ethyl','hyp_cata_ndi_ethyl')


def finite(value):return type(value) in (int,float) and math.isfinite(value)


def core_topology_plate(rows,out):
    import cairosvg
    by_id={r['id']:r for r in rows};cells=[]
    for gid,label in zip(IDS,('Conventional NDI / ethyl','Cata-NDI / ethyl')):
        row=by_id[gid];smiles=row['exact_grafted_smiles'];mol=Chem.MolFromSmiles(smiles)
        cells.append(dict(id=gid,label=label,smiles=smiles,formula=rdMolDescriptors.CalcMolFormula(mol),E1_vs_Fc_V=row.get('E1_vs_Fc_V') if row['potential_display_status']=='Admitted calculation' else None,E2_vs_Fc_V=row.get('E2_vs_Fc_V') if row['potential_display_status']=='Admitted calculation' else None,capacity_nominal_mAh_g=row['capacity_nominal_mAh_g'],potential_display_status=row['potential_display_status']))
    same_formula=cells[0]['formula']==cells[1]['formula']
    admitted=all(c['potential_display_status']=='Admitted calculation' and all(finite(c.get(k)) for k in ('E1_vs_Fc_V','E2_vs_Fc_V')) for c in cells)
    delta={k:cells[1][k]-cells[0][k] for k in ('E1_vs_Fc_V','E2_vs_Fc_V')} if admitted else {}
    parts=['<svg xmlns="http://www.w3.org/2000/svg" width="7.2in" height="4.1in" viewBox="0 0 1080 615">','<rect width="1080" height="615" fill="white"/>','<g font-family="DejaVu Sans,Arial,sans-serif" fill="#202a35">']
    for i,cell in enumerate(cells):
        x=540*i;center=x+270
        parts.append(f'<text x="{center}" y="38" text-anchor="middle" font-size="26" font-weight="bold">{html.escape(cell["label"])}</text>')
        parts.append(f'<g transform="translate({x+27},55) scale(1.35)">{molecule_svg(cell["smiles"])}</g>')
        potentials=' | '.join(f'E{j}: {cell[key]:.3f} V' if finite(cell[key]) else f'E{j}: pending' for j,key in enumerate(('E1_vs_Fc_V','E2_vs_Fc_V'),1))
        parts.append(f'<text x="{center}" y="441" text-anchor="middle" font-size="24">{potentials}</text>')
        parts.append(f'<text x="{center}" y="482" text-anchor="middle" font-size="23">Nominal capacity: {cell["capacity_nominal_mAh_g"]:.1f} mAh g⁻¹</text>')
    insight=(f'ΔE₁ = {delta["E1_vs_Fc_V"]:+.3f} V   ·   ΔE₂ = {delta["E2_vs_Fc_V"]:+.3f} V' if admitted else 'Potential comparison pending source-current admission')
    parts.extend(['<rect x="20" y="515" width="1040" height="55" rx="6" fill="#eef5f5"/>',f'<text x="540" y="551" text-anchor="middle" font-size="26" fill="#0e7c86">{insight}</text>',f'<text x="540" y="602" text-anchor="middle" font-size="19">{"Same molecular formula. " if same_formula else ""}Potentials vs Fc/Fc⁺ in MeCN. Local graft models; polymer effects are separate.</text>','</g></svg>'])
    svg='\n'.join(parts);(out/'core_topology_comparison.svg').write_text(svg)
    cairosvg.svg2pdf(bytestring=svg.encode(),write_to=str(out/'core_topology_comparison.pdf'))
    cairosvg.svg2png(bytestring=svg.encode(),write_to=str(out/'core_topology_comparison.png'),dpi=600)
    from PIL import Image
    with Image.open(out/'core_topology_comparison.png') as raster:raster.save(out/'core_topology_comparison.png',dpi=(600,600))
    report=dict(status='source_bound_core_topology_comparison',structures=cells,same_formula=same_formula,comparison_published=admitted,delta=delta,scope='computed capped local models; no vinyl or polymer properties transferred')
    (out/'core_topology_comparison.json').write_text(json.dumps(report,indent=2)+'\n')
    return report
