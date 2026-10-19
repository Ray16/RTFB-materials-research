#!/usr/bin/env python
"""Update the same all-candidate figures and comparison browser from the verified feed.

Default outputs: results/candidates/figures/. No scientific jobs are submitted.
Plotly is optional; matplotlib supplies the figures when it is unavailable.
"""
import argparse
import csv
import fcntl
import hashlib
import html
import json
import math
import os
from pathlib import Path
import sys
import tempfile
from functools import lru_cache

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(Path(__file__).resolve().parent))
from verified_overlay import load

ENERGY_FIELDS = ('E1_vs_Fc_V', 'E2_vs_Fc_V', 'wave_spacing_V',
                 'dG_disproportionation_kJ_mol', 'accessible_capacity_mAh_g',
                 'lambda_i_eV', 'lambda_outer_eV')
COLORS = ['#0072B2', '#D55E00', '#009E73', '#CC79A7', '#E69F00', '#56B4E9']
FAMILY_STYLE = {
    'imide (n-type)': ('Diimides', '#0072B2', 'o'),
    'quinone (n-type)': ('Quinones', '#D55E00', 's'),
    'pyridine-multi-e': ('Viologens / pyridines', '#009E73', '^'),
    'mono-imide (n-type)': ('Monoimides', '#CC79A7', 'D'),
    'phenothiazine (p-type, two-electron proposal)': ('Phenothiazines', '#E69F00', 'v'),
    'amine (p-type)': ('Amine catholytes', '#56B4E9', 'P'),
}


def comparison_rows(payload, rows, admitted):
    """Keep every candidate; source-invalid numeric properties cannot become points."""
    original = {r['id'] for r in payload['original39']}
    valid = {r['id'] for r in admitted}
    if len({r['id'] for r in rows}) != len(rows):
        raise ValueError('Candidate IDs must be unique across cohorts')
    result = []
    for source in rows:
        row = dict(source, cohort='Original' if source['id'] in original else 'Additional')
        if row['id'] not in valid:
            for key in ENERGY_FIELDS:
                row[key] = None
        if row.get('model_primary_admissible') is False:
            row['accessible_capacity_mAh_g'] = None
        row['potential_display_status'] = ('Admitted calculation' if row['id'] in valid
                                           else 'Pending / excluded / source changed')
        result.append(row)
    return sorted(result, key=lambda r: (r['family'], r['id']))


def finite(value):
    return type(value) in (int, float) and math.isfinite(value)


@lru_cache(maxsize=1)
def label_module():
    from importlib.util import module_from_spec, spec_from_file_location
    # Use the same labels as the slide tables, including exact aryl attachment sites.
    spec = spec_from_file_location('comparison_slide_labels', ROOT/'scripts/pipeline/update_verified_slides.py')
    module = module_from_spec(spec); spec.loader.exec_module(module)
    return module


def short_name(row, labels):
    return label_module().candidate_label(row, labels)


def figures(rows, out, stamp, report):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.backends.backend_pdf import PdfPages
    families = sorted({r['family'] for r in rows})
    styles = {name: FAMILY_STYLE.get(name, (name, '#555555', 'X')) for name in families}
    colors = {name: styles[name][1] for name in families}
    labels_path = ROOT/'results/today_2026-10-11_verified_collection/candidate_display_names.json'
    labels = json.loads(labels_path.read_text()) if labels_path.exists() else {}
    names = {r['id']: short_name(r, labels) for r in rows}
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 8,
        'axes.labelsize': 9, 'axes.titlesize': 9, 'xtick.labelsize': 7,
        'ytick.labelsize': 7, 'pdf.fonttype': 42, 'ps.fonttype': 42,
        'svg.fonttype': 'none', 'axes.spines.top': False, 'axes.spines.right': False,
        'savefig.bbox': None, 'font.weight': 'normal', 'axes.titleweight': 'normal',
        'axes.labelweight': 'normal'})
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 4.5), sharey=True)
    fig.subplots_adjust(left=.1, right=.98, top=.9, bottom=.26, wspace=.15)
    plotted = {}
    for ax, field, title in zip(axes, ('E1_vs_Fc_V', 'E2_vs_Fc_V'), ('First event', 'Second event')):
        plotted[field] = []
        for family in families:
            subset = [r for r in rows if r['family'] == family and finite(r.get(field))
                      and finite(r.get('capacity_nominal_mAh_g'))]
            plotted[field].extend(r['id'] for r in subset)
            ax.scatter([r[field] for r in subset], [r['capacity_nominal_mAh_g'] for r in subset],
                       facecolors=['none' if r.get('model_primary_admissible') is False else colors[family] for r in subset],
                       marker=styles[family][2], s=28, alpha=.85,
                       label=styles[family][0], edgecolors=colors[family], linewidths=.7)
        ax.set(xlabel='Potential (V vs Fc/Fc⁺)', title=title)
        ax.grid(alpha=.18)
        for row in rows:
            if row['id'] not in report['highlighted_ids']:continue
            ax.scatter(row[field],row['capacity_nominal_mAh_g'],marker='*',s=85,color='#111',zorder=5)
            ax.annotate(names[row['id']],(row[field],row['capacity_nominal_mAh_g']),
                        xytext=(5,5),textcoords='offset points',fontsize=6.5)
    axes[0].set_ylabel('Nominal capacity (mAh g⁻¹)')
    for letter, ax in zip('ab', axes):ax.text(-.12, 1.08, letter, transform=ax.transAxes, fontweight='bold', fontsize=11)
    handles = [Line2D([], [], color=colors[f], marker=styles[f][2], linestyle='', label=styles[f][0]) for f in families]
    handles.append(Line2D([],[],color='#444',marker='o',markerfacecolor='none',linestyle='',label='Open: model selection hold'))
    handles.append(Line2D([],[],color='#111',marker='*',linestyle='',label='Selected anolyte descriptor leaders'))
    fig.legend(handles=handles, loc='lower center', bbox_to_anchor=(.5,.06), ncol=3, fontsize=7.5, frameon=False)
    for suffix in ('png', 'pdf', 'svg'):fig.savefig(out/f'redox_capacity.{suffix}', dpi=600)
    plt.close(fig)

    fig, (left, right) = plt.subplots(1, 2, figsize=(14, max(8, .245*len(rows))),
                                    sharey=True, gridspec_kw={'width_ratios': [3, 1]})
    for i, row in enumerate(rows):
        values = [row.get(k) for k in ('E1_vs_Fc_V', 'E2_vs_Fc_V')]
        for j, value in enumerate(values):
            if finite(value):left.scatter(value, i, marker='o' if j == 0 else 's', color=colors[row['family']], s=22)
        if all(finite(v) for v in values):left.plot(values, [i, i], color=colors[row['family']], alpha=.5)
        if not any(finite(v) for v in values):left.text(.02, i, 'pending / excluded', transform=left.get_yaxis_transform(), color='#777', fontsize=7)
        if finite(row.get('capacity_nominal_mAh_g')):
            right.barh(i, row['capacity_nominal_mAh_g'], color=colors[row['family']], alpha=.75)
    left.set_yticks(range(len(rows)), [names[r['id']] for r in rows], fontsize=7)
    left.invert_yaxis();left.set_xlabel('Potential (V vs Fc/Fc⁺, MeCN) • circle E₁, square E₂')
    right.set_xlabel('Nominal capacity\n(mAh/g)')
    for ax in (left, right):ax.grid(axis='x', alpha=.18)
    fig.suptitle(f'All-candidate comparison • {len(rows)} exact structures\n{stamp} • capacity excludes electrolyte and polymer', fontsize=12)
    fig.subplots_adjust(left=.35, right=.98, top=.95, bottom=.04, wspace=.09)
    fig.savefig(out/'redox_ladders.png', dpi=300)
    plt.close(fig)

    # Readable supplementary pages instead of shrinking 75 labels into a paper panel.
    with PdfPages(out/'redox_ladders.pdf') as document:
        for offset in range(0, len(rows), 18):
            chunk = rows[offset:offset+18]
            fig, (ax, cap) = plt.subplots(1, 2, figsize=(7.2, 6.5), sharey=True,
                gridspec_kw={'width_ratios':[3,1]})
            for i, row in enumerate(chunk):
                values=[row.get(k) for k in ('E1_vs_Fc_V','E2_vs_Fc_V')]
                for j, value in enumerate(values):
                    if finite(value):ax.scatter(value,i,marker='o' if j==0 else 's',color=colors[row['family']],s=22)
                if all(finite(v) for v in values):ax.plot(values,[i,i],color=colors[row['family']],alpha=.6)
                if not any(finite(v) for v in values):ax.text(.02,i,'pending / excluded',transform=ax.get_yaxis_transform(),fontsize=7,color='#777')
                if finite(row.get('capacity_nominal_mAh_g')):cap.barh(i,row['capacity_nominal_mAh_g'],color=colors[row['family']],alpha=.8)
            ax.set_yticks(range(len(chunk)),[names[r['id']] for r in chunk],fontsize=7)
            ax.invert_yaxis();ax.set_xlabel('Potential (V vs Fc/Fc⁺)');cap.set_xlabel('Nominal capacity\n(mAh g⁻¹)')
            fig.suptitle(f'All-candidate comparison · {offset+1}–{offset+len(chunk)} of {len(rows)}',fontsize=10)
            fig.subplots_adjust(left=.4,right=.98,top=.92,bottom=.1,wspace=.15)
            document.savefig(fig);plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.2, 3.8), constrained_layout=True)
    complete = [sum(r['family'] == f and r['potential_display_status'] == 'Admitted calculation' for r in rows) for f in families]
    pending = [sum(r['family'] == f for r in rows)-n for f, n in zip(families, complete)]
    family_labels=[styles[f][0] for f in families]
    ax.barh(family_labels, complete, label='Admitted calculation', color='#0072B2')
    ax.barh(family_labels, pending, left=complete, label='Pending / excluded', color='#c9cdd2',hatch='///',edgecolor='#999',linewidth=.3)
    for i, (n, p) in enumerate(zip(complete, pending)):ax.text(n+p+.2, i, f'{n}/{n+p}')
    ax.set(xlabel='Candidates', title=f'Potential calculation coverage · {sum(complete)}/{len(rows)} candidates')
    ax.set_xlim(0, max(a+b for a, b in zip(complete, pending))*1.15)
    ax.legend(loc='lower right',frameon=False,fontsize=8);ax.grid(axis='x', alpha=.15)
    for suffix in ('png', 'pdf', 'svg'):fig.savefig(out/f'coverage.{suffix}', dpi=600)
    plt.close(fig)
    return colors, names, plotted


def browser(rows, colors, names, stamp, report):
    from rdkit import Chem
    from rdkit.Chem.Draw import rdMolDraw2D
    structures = {}
    for row in rows:
        mol = Chem.MolFromSmiles(row['exact_grafted_smiles'])
        if mol is None:raise ValueError('Invalid exact structure: '+row['id'])
        drawer = rdMolDraw2D.MolDraw2DSVG(420, 230)
        drawer.DrawMolecule(mol);drawer.FinishDrawing()
        structures[row['id']] = drawer.GetDrawingText().split('?>')[-1]
    try:
        import plotly.graph_objects as go
        from plotly.subplots import make_subplots
        fig = make_subplots(rows=1, cols=2, subplot_titles=['First event', 'Second event'])
        for column, field in enumerate(('E1_vs_Fc_V', 'E2_vs_Fc_V'), 1):
            for family in sorted(colors):
                subset = [r for r in rows if r['family'] == family and finite(r.get(field))]
                fig.add_trace(go.Scatter(x=[r[field] for r in subset],
                    y=[r['capacity_nominal_mAh_g'] for r in subset], mode='markers',
                    marker=dict(color=colors[family],
                        size=[14 if r['id'] in report['highlighted_ids'] else 9 for r in subset],
                        symbol=['star' if r['id'] in report['highlighted_ids'] else 'circle-open' if r.get('model_primary_admissible') is False else 'circle' for r in subset]), name=family,
                    legendgroup=family, showlegend=column == 1,
                    customdata=[[names[r['id']], r['id']] for r in subset],
                    hovertemplate='%{customdata[0]}<br>%{customdata[1]}<br>%{x:.3f} V<br>%{y:.1f} mAh/g<extra></extra>'), row=1, col=column)
            fig.update_xaxes(title_text='V vs Fc/Fc⁺ in MeCN', row=1, col=column)
            fig.update_yaxes(title_text='Nominal capacity (mAh/g)', row=1, col=column)
        fig.update_layout(template='plotly_white', height=510, margin=dict(l=60, r=25, t=45, b=90),
                          legend=dict(orientation='h', y=-.25), font=dict(family='system-ui'))
        plot = fig.to_html(full_html=False, include_plotlyjs=True, div_id='redox-chart', config=dict(responsive=True))
        interactive = True
    except ImportError:
        plot = '<img src="redox_capacity.png" alt="Redox potentials versus nominal capacity" style="width:100%">'
        interactive = False
    table = []
    for row in rows:
        display = dict(row, short_name=names[row['id']], structure_svg=structures[row['id']])
        table.append(display)
    data = json.dumps(table, allow_nan=False).replace('</', '<\\/')
    page = (ROOT/'scripts/plotting/candidates/comparison_browser.html').read_text()
    page = page.replace('__STAMP__', html.escape(stamp)).replace('__PLOT__', plot).replace('__DATA__', data)
    return page, interactive


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--feed', type=Path, default=ROOT/'results/candidates/current.json')
    parser.add_argument('--out', type=Path, default=ROOT/'results/candidates/figures')
    args = parser.parse_args(argv);args.out.mkdir(parents=True, exist_ok=True)
    lock = (args.out/'.update.lock').open('a+')
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    producers=[Path(__file__),Path(__file__).with_name('candidate_highlights.py'),
        Path(__file__).with_name('comparison_browser.html'),Path(__file__).with_name('verified_overlay.py'),
        ROOT/'scripts/pipeline/update_verified_slides.py']
    producer_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in producers}
    payload, raw, admitted = load(args.feed)
    if not payload or not raw:raise RuntimeError('No source-current audited candidate snapshot; collection preserved')
    _,_,comparison_admitted=load(args.feed,include_controls=True)
    rows = comparison_rows(payload, raw, admitted)
    from candidate_highlights import highlights,matched_highlights,focused_figures
    report=highlights(rows);matches=matched_highlights(payload,comparison_admitted)
    stamp = payload['as_of_Chicago'];snapshot = payload['snapshot_path']
    with tempfile.TemporaryDirectory(prefix='.build-', dir=args.out) as temp:
        stage = Path(temp)
        colors, names, points = figures(rows, stage, stamp, report)
        for row in comparison_admitted:
            if row['id'] not in names:names[row['id']]=short_name(row,{})
        focused_figures(rows,report,matches,stage,names,colors)
        (stage/'highlights.json').write_text(json.dumps(dict(report,matched_contrasts=matches),indent=2)+'\n')
        slide_frames=[]
        for title,filename,caption in [
            ('Candidate descriptor leaders: potential and capacity','redox_capacity.png',
             'Stars select anolyte descriptor leaders on a partial potential--capacity front. Kinetics and polymer validation remain pending.'),
            ('Synthesis accessibility is a separate consideration','synthesis_tradeoff.png',
             'SA is a heuristic annotation. The computed methylbenzyl caps lack the vinyl group required for monomer-first polymerization.'),
            ('New candidate contrasts against matched controls','new_candidate_shifts.png',
             'One contrast per core with the strongest potential decrease shared by both events. A shift does not establish polymer performance.')]:
            slide_frames.extend([r'\begin{frame}{'+title+'}',r'\centering',
                r'\includegraphics[width=.97\linewidth,height=.72\textheight,keepaspectratio]{../results/candidates/figures/'+filename+'}',
                r'\par\smallskip\raggedright\scriptsize '+caption,
                r'\vfill\tiny Audited model snapshot: '+stamp.replace('_',r'\_')+'.',r'\end{frame}'])
        (stage/'candidate_comparison_collection.tex').write_text('\n'.join(slide_frames)+'\n')
        page, interactive = browser(rows, colors, names, stamp, report)
        (stage/'index.html').write_text(page)
        (stage/'all_candidates.json').write_text(json.dumps(rows, indent=2, allow_nan=False)+'\n')
        (stage/'captions.md').write_text(
            '# Figure captions\n\n'
            f'Audited snapshot: {stamp}. {len(rows)} registered exact grafts, '
            f'{len(admitted)} with admitted potential calculations.\n\n'
            '**Potential–capacity map.** First (a) and second (b) declared redox-event '
            'potentials versus Fc/Fc+ in MeCN and nominal gravimetric capacity of the '
            'p-methylbenzyl-capped local graft model. Each point is a candidate; colors and shapes identify '
            'families. Open symbols mark model-selection holds on qualified capacity. '
            'Nominal capacity follows the declared electron count and molecular mass, '
            'excluding counterions, electrolyte, solvent and polymer. Missing potential '
            'endpoints are omitted, not set to zero. Anolyte and catholyte roles must '
            'be assessed separately; no overall ranking or cell voltage follows.\n\n'
            'Stars highlight the two highest nominal capacities and most negative first '
            'reduction on the partial anolyte potential–capacity front. This exact '
            'descriptor front omits kinetics and candidate-specific uncertainty; it '
            'does not replace the qualified primary Pareto policy.\n\n'
            '**Synthesis tradeoff.** Nominal molecular capacity versus heuristic SA. '
            'Colored points lie on the partial potential–capacity front; stars identify '
            'the same selected descriptor leaders. Gray points are other or pending '
            'designs. SA is an annotation, not a measured synthesis outcome or ranking '
            'objective. The model caps do not contain polymerizable vinyl groups.\n\n'
            '**New candidate shifts.** One matched contrast per core with the largest '
            'potential decrease shared by both redox events, among source-current '
            'two-electron anolyte calculations passing the available window/integrity '
            'bookkeeping. Both potentials are differences of exact local-model '
            'calculations under the same protocol. A favorable shift does not '
            'establish polymer synthesis, stability, loading or cycling performance.\n\n'
            '**All-candidate ladders (supplementary comparison).** Circles and squares '
            'show first and second declared potential events. Lines connect endpoints '
            'for the same structure and do not represent measured reaction paths. '
            'Adjacent bars show nominal molecular capacity. Pending or excluded '
            'endpoints remain blank. The PDF uses readable supplementary pages; the '
            'PNG is an overview of the complete collection.\n\n'
            '**Family coverage.** Candidate counts with source-current admitted '
            'potential calculations (blue) and pending/excluded calculations (hatched '
            'gray). Admission applies to the potential calculation, not to experimental '
            'accuracy, kinetics, solubility or durability.\n\n'
            'All potentials use r2SCAN-D4/SVP(D)/SMD geometries, wB97M-V/TZVPD/SMD '
            'active energies and independently checked native internal-mode '
            'thermochemistry at 298.15 K with the matched Fc reference. These are '
            'model predictions for isolated exact grafts. No experimental error bars '
            'are assigned: general exact-graft error bounds are unestablished. The '
            'phenothiazine second-wave model discrepancy remains unresolved. '
            'Wave spacing and disproportionation derive from the same endpoint '
            'energies and are not independent evidence. Source-parent measurements '
            'are not transferred to grafts. Solution reorganization values remain '
            'blank when unqualified.\n')
        fields = list(rows[0])
        with (stage/'all_candidates.csv').open('w', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=fields);writer.writeheader();writer.writerows(rows)
        current, raw_after, admitted_after = load(args.feed)
        _,_,comparison_after=load(args.feed,include_controls=True)
        if (not current or current['snapshot_path'] != snapshot
                or comparison_rows(current, raw_after, admitted_after) != rows
                or comparison_after!=comparison_admitted
                or any(hashlib.sha256(Path(p).read_bytes()).hexdigest()!=h for p,h in producer_sha256.items())):
            raise RuntimeError('Source changed during rendering; collection preserved, rerun update')
        manifest = dict(status='updated_verified_comparison_collection', as_of_Chicago=stamp,
            snapshot=snapshot, candidates=len(rows), original=sum(r['cohort']=='Original' for r in rows),
            additional=sum(r['cohort']=='Additional' for r in rows),
            admitted_potential_candidates=len(admitted), plotted_candidate_ids=points,
            interactive_plotly=interactive, structure_views=len(rows),
            producer_sha256=producer_sha256,
            publication_style=dict(width_inches=7.2,main_PNG_dpi=600,vector_formats=['pdf','svg'],
                embedded_TrueType_fonts=True,grayscale_redundancy='family marker shapes and pending hatching',
                supplementary_rows_per_PDF_page=18),
            exact_graft_experimental_accuracy_validated=False,
            files={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in stage.iterdir()})
        (stage/'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
        # Publish the manifest last so it always identifies the completed collection.
        for path in sorted(stage.iterdir(), key=lambda p:p.name=='manifest.json'):
            os.replace(path, args.out/path.name)
        slide_source=ROOT/'slides/generated/candidate_comparison_collection.tex'
        slide_source.parent.mkdir(parents=True,exist_ok=True)
        temporary=slide_source.with_name('.candidate_comparison_collection.tmp')
        temporary.write_bytes((args.out/'candidate_comparison_collection.tex').read_bytes())
        temporary.replace(slide_source)
    print(json.dumps({k:v for k,v in manifest.items() if k not in ('files','plotted_candidate_ids')}))


if __name__ == '__main__':main()
