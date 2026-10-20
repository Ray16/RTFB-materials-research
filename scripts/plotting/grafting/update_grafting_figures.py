#!/usr/bin/env python
"""Maintain focused grafting figures and a searchable exact-structure browser."""
import hashlib
import html
import io
import json
import os
from pathlib import Path
import subprocess
import textwrap
import xml.etree.ElementTree as ET

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import rdChemReactions
from rdkit.Chem.Draw import rdMolDraw2D

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT/'results/grafting'
FIG = OUT/'figures'
COLORS = {'NH_OH_substitution':'#0072B2', 'N_quaternization':'#D55E00'}
LABELS = {'NH_OH_substitution':'N–H / O–H substitution (HCl reference)',
          'N_quaternization':'N quaternization (chloride reference)'}


def save(fig, name):
    for suffix in ('pdf', 'svg', 'png'):
        fig.savefig(FIG/f'{name}.{suffix}', dpi=350, facecolor='white')
    plt.close(fig)


def structure(smiles, destination):
    if destination.exists(): return
    mol = Chem.MolFromSmiles(smiles)
    drawer = rdMolDraw2D.MolDraw2DSVG(360, 210)
    drawer.drawOptions().bondLineWidth = 1.5
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    destination.write_text(drawer.GetDrawingText())


def main():
    FIG.mkdir(parents=True, exist_ok=True)
    structures = OUT/'structures'; structures.mkdir(exist_ok=True)
    feed = json.loads((OUT/'current.json').read_text())
    rows, status = feed['candidates'], feed['status']
    inventory = json.loads((OUT/'reactions.json').read_text())['reactions']
    assessment_path = OUT/'selection_assessment.json'
    assessment = json.loads(assessment_path.read_text()) if assessment_path.exists() else {}
    synthesis = assessment.get('candidates', {})
    plt.rcParams.update({'font.family':'DejaVu Sans', 'font.size':9, 'axes.labelsize':10,
                         'axes.spines.top':False, 'axes.spines.right':False,
                         'svg.fonttype':'none', 'pdf.fonttype':42})
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 5.6), layout='constrained')
    for ax, kind in zip(axes, COLORS):
        all_class = [r for r in rows if r['reaction_class']==kind]
        ready = sorted((r for r in all_class if r['deltaG_reference_kJ_mol'] is not None),
                       key=lambda r:r['deltaG_reference_kJ_mol'])[:12]
        nready = sum(r['deltaG_reference_kJ_mol'] is not None for r in all_class)
        ax.set_title(LABELS[kind]+f'\n{nready}/{len(all_class)} qualified', fontsize=10, pad=12)
        if ready:
            ys = list(range(len(ready)))
            ax.scatter([r['deltaG_reference_kJ_mol'] for r in ready], ys,
                       color=COLORS[kind], s=38, zorder=3)
            ax.set_yticks(ys, [r['id'].removeprefix('exp_').removeprefix('design_') for r in ready], fontsize=8)
            ax.invert_yaxis(); ax.axvline(0, color='.65', lw=.8)
            ax.grid(axis='x', alpha=.15); ax.set_axisbelow(True)
        else:
            ax.text(.5, .5, 'Awaiting qualified reaction components', ha='center', va='center',
                    transform=ax.transAxes, fontsize=10, color='.4', wrap=True)
            ax.set_yticks([]); ax.set_xticks([])
        ax.set_xlabel('Reference grafting ΔG (kJ mol⁻¹)')
    fig.suptitle('Local grafting thermodynamics · MeCN, 298.15 K', fontsize=13)
    fig.supxlabel('Up to 12 lowest reference ΔG values per reaction class. Base capture, barriers and polymer effects are separate.',
                  fontsize=8)
    save(fig, 'deltaG_ranking')

    fig, ax = plt.subplots(figsize=(7.2, 3.0), layout='constrained')
    groups = ['Qualified reactions','Pending components','Held by component QC']
    counts = [sum(r['status']=='qualified_reference_reaction' for r in rows),
              sum(r['status']=='pending_components' for r in rows),
              sum(r['status']=='held_failed_component_QC' for r in rows)]
    bars = ax.barh(groups, counts, color=['#009E73','#999999','#D55E00'], height=.5)
    ax.bar_label(bars, padding=5); ax.set_xlim(0, max(80, max(counts)+8))
    ax.invert_yaxis(); ax.set_xlabel('Candidates'); ax.set_title('Grafting calculation coverage')
    ax.spines['left'].set_visible(False); ax.tick_params(axis='y', length=0)
    fig.supxlabel(f"{status['qualified_components']}/{status['expected_new_components']} new components qualified; "
                  f"{status['preserved_component_failures']} component failures preserved", fontsize=9)
    save(fig, 'coverage')

    # Chemical schemes are static until exact reaction identities change.
    signature = 'vector_v2:'+hashlib.sha256((OUT/'reactions.json').read_bytes()).hexdigest()
    marker = FIG/'reaction_schemes.source'
    if not marker.exists() or marker.read_text()!=signature:
        chosen = ['exp_bq_oxy_methoxy', 'design_ndi_ethyl', 'exp_viologen_methoxyethyl']
        schemes = [next(r for r in inventory if r['id']==key) for key in chosen]
        svg = ['<svg xmlns="http://www.w3.org/2000/svg" width="1056" height="700" viewBox="0 0 1056 700">',
               '<rect width="1056" height="700" fill="white"/>',
               '<text x="24" y="32" font-family="DejaVu Sans" font-size="22">Local attachment reference reactions</text>']
        labels = ['Methoxy-benzoquinone · O attachment', 'Naphthalene diimide · N attachment',
                  'Viologen · N quaternization']
        for i, (r, label) in enumerate(zip(schemes, labels)):
            scheme = r['site_smiles']+'.'+r['precursor_smiles']+'>>'+r['product_smiles']+'.'+r['byproduct_smiles']
            rxn = rdChemReactions.ReactionFromSmarts(scheme, useSmiles=True)
            drawer = rdMolDraw2D.MolDraw2DSVG(1800, 300)
            drawer.drawOptions().bondLineWidth = 2.5
            drawer.DrawReaction(rxn); drawer.FinishDrawing()
            nested = ET.fromstring(drawer.GetDrawingText())
            nested.attrib.update(x='0', y=str(70+i*195), width='1056', height='176')
            svg.append(f'<text x="24" y="{68+i*195}" font-family="DejaVu Sans" font-size="15">{html.escape(label)}</text>')
            svg.append(ET.tostring(nested, encoding='unicode'))
        svg.append('<text x="24" y="672" font-family="DejaVu Sans" font-size="12">MeCN, 298.15 K · capped-site models · base capture and polymer environment assessed separately</text></svg>')
        (FIG/'reaction_schemes.svg').write_text(''.join(svg))
        # CairoSVG lives in the existing isolated figure environment. Both
        # molecular paths and text remain vector objects in the PDF.
        renderer = ROOT/'.venv-figures/bin/python'
        code = ('import cairosvg,sys;'
                'cairosvg.svg2pdf(url=sys.argv[1],write_to=sys.argv[2]);'
                'cairosvg.svg2png(url=sys.argv[1],write_to=sys.argv[3],output_width=3850)')
        subprocess.run([str(renderer), '-c', code, str(FIG/'reaction_schemes.svg'),
                        str(FIG/'reaction_schemes.pdf'), str(FIG/'reaction_schemes.png')], check=True)
        marker.write_text(signature)

    cards = []
    for r in rows:
        product = structures/(r['id']+'.svg')
        structure(r['product_smiles'], product)
        value = r['deltaG_reference_kJ_mol']
        value_text = f'{value:+.1f} kJ mol⁻¹' if value is not None else 'ΔG pending'
        decision = synthesis.get(r['id'], {})
        selection_text = ('<p class="hold"><strong>'+html.escape(decision.get('display_label', decision.get('status', 'Route review required')))
                          +'</strong><br>'+html.escape(decision.get('reason', ''))+'</p>') if decision else ''
        cards.append('<article data-search="'+html.escape(r['id']+' '+r['name']+' '+r['family']+' '+r['status']+' '+decision.get('status','unassessed'), quote=True)+'">'
            '<h3>'+html.escape(r['id'])+'</h3><img loading="lazy" alt="Exact graft-site product structure" src="structures/'+html.escape(product.name)+'">'
            '<p>'+html.escape(r['name'])+'</p><strong>'+value_text+'</strong><p class="muted">'+html.escape(r['status'])+'</p>'+selection_text+
            '<details><summary>Reaction and QC</summary><p>'+html.escape(LABELS[r['reaction_class']])+'</p>'
            '<p>Missing: '+html.escape(r['missing_components'] or 'none')+'</p><p>Held: '+html.escape(r['failed_components'] or 'none')+'</p>'
            '<p>Product records: <code>'+html.escape(r['product_state_path'])+'</code></p></details></article>')
    page = '''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Grafting free energies</title><style>
body{font:16px/1.5 system-ui,sans-serif;color:#203040;margin:30px auto;max-width:1250px;padding:0 20px;background:#f8fafb}
h1{font-size:30px}a{color:#0072B2}.muted{color:#627080}input{padding:12px;width:min(95%,600px);font:inherit;border:1px solid #bbc8d0;border-radius:6px}
.hold{background:#fff4df;border-left:4px solid #b87500;padding:10px}
.plots{display:grid;grid-template-columns:2fr 1fr;gap:20px;align-items:center}.plots img,.scheme{width:100%}
.cards{display:grid;grid-template-columns:repeat(auto-fit,minmax(285px,1fr));gap:18px;margin-top:24px}
article{background:white;border:1px solid #dbe2e8;border-radius:9px;padding:18px;overflow-wrap:anywhere}article img{width:100%;height:185px;object-fit:contain}
h3{font-size:16px}details{font-size:13px}code{font-size:11px}@media(max-width:750px){.plots{display:block}}
</style><h1>Grafting free energies</h1><p><a href="../README.md">All results</a> · <a href="README.md">Methods and scope</a> · <a href="current.csv">Download table</a> · <a href="status.json">Status and provenance</a> · <a href="base_assisted/index.html">Base-assisted pilot</a></p>
'''
    page += f"<p><strong>{status['completed_grafting_deltaG']}/75 reference reaction ΔG values qualified.</strong> "
    page += f"{status['qualified_components']}/{status['expected_new_components']} new components qualified. Reload this page to see collected updates.</p>"
    page += '<p>Acetonitrile, 298.15 K. N–H/O–H reactions retain HCl explicitly; N quaternization retains chloride. '
    page += 'Compare the two reaction classes separately. The completed <a href="base_assisted/index.html">DBU/DIPEA pilot</a> reports base capture separately. These local-site predictions do not establish grafting rates, yields or polymer performance.</p>'
    if synthesis:
        page += '<p class="hold"><strong>Base required; BQ candidates retained for evaluation.</strong> The modeled DBU/DIPEA routes are favorable for both BQ designs and AQ-1-O. '
        page += 'AQ-2-O/6-OMe remains unresolved because of precursor QC. Positive HCl-reference values alone do not exclude a candidate. '
        page += 'Attachment selectivity, redox retention and yield still need review. <a href="selection_assessment.json">Selection assessment</a></p>'
    page += '<div class="plots"><a href="figures/deltaG_ranking.pdf"><img src="figures/deltaG_ranking.png" alt="Reference grafting free energies; pending until qualified"></a>'
    page += '<a href="figures/coverage.pdf"><img src="figures/coverage.png" alt="Qualified, pending and QC-held reaction counts"></a></div>'
    page += '<a href="figures/reaction_schemes.pdf"><img class="scheme" src="figures/reaction_schemes.png" alt="Three balanced attachment reaction schemes"></a>'
    page += '<h2>All candidates and exact molecular structures</h2><input id="search" type="search" placeholder="Filter by candidate, family or status" aria-label="Filter candidates">'
    page += '<div class="cards">'+''.join(cards)+'</div><script>document.getElementById("search").addEventListener("input",function(){const q=this.value.toLowerCase();document.querySelectorAll("article").forEach(c=>{c.hidden=!c.dataset.search.toLowerCase().includes(q)})})</script></html>'
    temporary = OUT/'index.html.tmp'; temporary.write_text(page); temporary.replace(OUT/'index.html')
    artifacts = sorted([*FIG.glob('*.pdf'), *FIG.glob('*.png'), *FIG.glob('*.svg'), *structures.glob('*.svg')])
    atomic = dict(source_current_sha256=hashlib.sha256((OUT/'current.json').read_bytes()).hexdigest(),
                  figure_producer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  synthesis_assessment_sha256=hashlib.sha256(assessment_path.read_bytes()).hexdigest() if assessment_path.exists() else None,
                  qualified_reactions=status['completed_grafting_deltaG'],
                  outputs={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in artifacts})
    (FIG/'manifest.json').write_text(json.dumps(atomic, indent=2)+'\n')
    print('Maintained grafting browser and three figure sets refreshed')


if __name__=='__main__': main()
