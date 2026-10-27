#!/usr/bin/env python3
"""Plot checked frozen local-motif electronic interactions versus separation."""
import argparse
import hashlib
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


EV_KJ_MOL=1.602176634e-19*6.02214076e23/1000


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''):h.update(block)
    return h.hexdigest()


def validated_points(path):
    panel=json.loads(path.read_text())
    if (panel['status']!='all_twelve_frozen_local_interaction_triplets_independently_reconstructed'
            or panel['completed_triplets']!=12 or panel['checked_electronic_geometries']!=36
            or any(panel[k] for k in ('free_energies_computed','electronic_coupling_computed',
                                      'conformer_convergence_established','basis_superposition_correction_computed'))):
        raise ValueError('Local electronic panel is incomplete or its scope differs')
    for source,digest in panel['source_sha256'].items():
        if sha(source)!=digest:raise ValueError('Source changed: '+source)
    grouped={}
    for row in panel['states']:
        if row['missing_states'] or row['free_energy_computed'] or row['basis_superposition_correction_computed']:
            raise ValueError('Incomplete or overstated local term')
        energies=row['exact_state_energies_eV'];actual=math.fsum([energies['AB']['acetonitrile'],
            -energies['A']['acetonitrile'],-energies['B']['acetonitrile']])*EV_KJ_MOL
        expected=row['interaction_energy_kJ_mol_by_medium']['acetonitrile']
        if not math.isfinite(actual) or abs(actual-expected)>1e-7:raise ValueError('Interaction does not reconstruct')
        grouped.setdefault(row['parent_key'],[]).append(row)
    expected={'occ_2_3__pose0','occ_2_3__pose1','occ_2_4__pose0','occ_2_4__pose1'}
    if set(grouped)!=expected:raise ValueError('Expected four sampled parent geometries')
    points=[]
    for parent,rows in sorted(grouped.items()):
        if sorted(round(r['cap_length_A'],2) for r in rows)!=[1.07,1.09,1.11]:raise ValueError('Incomplete cap-length comparison')
        distance=rows[0]['parent_contact']['centroid_distance_A']
        if not math.isfinite(distance) or any(abs(r['parent_contact']['centroid_distance_A']-distance)>1e-9 for r in rows):
            raise ValueError('Parent geometry differs between cap lengths')
        central=next(r for r in rows if round(r['cap_length_A'],2)==1.09)
        values=[r['interaction_energy_kJ_mol_by_medium']['acetonitrile'] for r in rows]
        points.append(dict(parent_key=parent,spacing='adjacent' if parent.startswith('occ_2_3') else 'one spacer',
            centroid_distance_A=distance,central_cap_A=1.09,interaction_kJ_mol=central['interaction_energy_kJ_mol_by_medium']['acetonitrile'],
            cap_range_kJ_mol=[min(values),max(values)],cap_lengths_A=[1.07,1.09,1.11]))
    return panel,points


def plot(panel_path,prefix):
    panel,points=validated_points(panel_path)
    outputs=[prefix.with_suffix(ext) for ext in ('.pdf','.png','.json')]
    if any(p.exists() for p in outputs):raise FileExistsError('Preserve existing figures and evidence')
    prefix.parent.mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({'font.size':15,'axes.labelsize':15,'xtick.labelsize':13,'ytick.labelsize':13,'legend.fontsize':13})
    fig,ax=plt.subplots(figsize=(6.5,3.65),layout='constrained')
    styles={'adjacent':('#0072B2','o'),'one spacer':('#D55E00','s')}
    for spacing,(color,marker) in styles.items():
        group=[p for p in points if p['spacing']==spacing]
        for index,point in enumerate(group):
            y=point['interaction_kJ_mol'];lo,hi=point['cap_range_kJ_mol']
            ax.errorbar(point['centroid_distance_A'],y,yerr=[[y-lo],[hi-y]],fmt=marker,color=color,
                        markersize=8,capsize=3,linewidth=1.1,label=spacing.title() if index==0 else None)
    ax.axhline(0,color='0.55',linewidth=.8,zorder=0)
    ax.set_xlim(0,23);ax.set_ylim(-31,5)
    ax.set_xlabel('Core-centroid separation (Å)')
    ax.set_ylabel(r'$E_{\mathrm{int}}$ (kJ/mol)')
    ax.spines[['top','right']].set_visible(False);ax.legend(frameon=False,loc='lower right')
    fig.savefig(outputs[0],metadata={'Creator':'redox workflow','CreationDate':None,'ModDate':None})
    fig.savefig(outputs[1],dpi=220,metadata={'Software':'redox workflow'});plt.close(fig)
    result=dict(status='checked_frozen_local_motif_contact_figure',points=points,source_panel=str(panel_path.resolve()),
        source_panel_sha256=sha(panel_path),source_files_verified=len(panel['source_sha256']),producer_sha256=sha(__file__),
        outputs={p.name:sha(p) for p in outputs[:2]},
        scope='Four sampled phthalimide capped-fragment poses in the MeCN model; error bars vary cap length only, not uncertainty or populations.',
        quantum_calculations_run=False,association_free_energy_inferred=False,electronic_coupling_inferred=False)
    with outputs[2].open('x') as stream:json.dump(result,stream,indent=2);stream.write('\n')
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--panel',type=Path,required=True);p.add_argument('--out-prefix',type=Path,required=True);a=p.parse_args()
    result=plot(a.panel.resolve(),a.out_prefix.resolve());print(json.dumps(dict(status=result['status'],points=len(result['points']))))
