#!/usr/bin/env python
"""Build the focused deck, recording compiler checks separately from visual review."""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import time

from pypdf import PdfReader

ROOT=Path(__file__).resolve().parents[2]


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def validate_tempo(tempo_evidence):
    tempo=json.loads(tempo_evidence.read_text())
    if (tempo['status']!='TEMPO_exact_support_route_pass_local_potential_checked_material_pending'
            or tempo['prospective_route_status']!='proceed' or tempo['all_axis_status']!='pending'
            or len(tempo['pending_axes'])!=9 or tempo['actual_route_deltaG_computed']
            or abs(tempo['modeled_oxidation_V_vs_Fc']['50']-.5787750407616841)>1e-9):
        raise ValueError('Unqualified TEMPO route or local potential')
    for name,expected in tempo['source_artifact_sha256'].items():
        if sha(tempo_evidence.parent/name)!=expected:raise ValueError('Changed TEMPO evidence: '+name)
    figure=ROOT/tempo['figure_pdf']
    if sha(figure)!=tempo['figure_pdf_sha256']:raise ValueError('Changed TEMPO figure')
    return tempo,figure


def validate_local_contacts():
    panel=ROOT/'results/reviews/2026-10-26/local_motif_complete/summary.json'
    metadata=ROOT/'results/figures/polymer_models/local_motif_contacts.json'
    figure=ROOT/'results/figures/polymer_models/local_motif_contacts.pdf'
    plotter=ROOT/'scripts/plotting/polymer_models/plot_local_motif_contacts.py'
    p=json.loads(panel.read_text());m=json.loads(metadata.read_text())
    if (p['status']!='all_twelve_frozen_local_interaction_triplets_independently_reconstructed'
            or p['completed_triplets']!=12 or p['checked_electronic_geometries']!=36
            or any(p[k] for k in ('free_energies_computed','electronic_coupling_computed',
                                  'conformer_convergence_established','basis_superposition_correction_computed'))
            or m['status']!='checked_frozen_local_motif_contact_figure' or len(m['points'])!=4
            or m['source_panel_sha256']!=sha(panel) or m['producer_sha256']!=sha(plotter)):
        raise ValueError('Local-contact figure or its scientific scope differs')
    for path,digest in p['source_sha256'].items():
        if sha(path)!=digest:raise ValueError('Changed local-contact source: '+path)
    pdf_hashes=[value for path,value in m['outputs'].items() if Path(path).name==figure.name]
    if len(pdf_hashes)!=1 or pdf_hashes[0]!=sha(figure):raise ValueError('Changed local-contact figure')
    rows={row['parent_key']:row for row in m['points']}
    if set(rows)!={'occ_2_3__pose0','occ_2_3__pose1','occ_2_4__pose0','occ_2_4__pose1'}:
        raise ValueError('Local-contact pose coverage differs')
    groups={key:[] for key in rows}
    for state in p['states']:
        if state['parent_key'] not in groups or state['missing_states']:
            raise ValueError('Local-contact source pose differs')
        energies=state['exact_state_energies_eV'];actual=math.fsum((energies['AB']['acetonitrile'],
            -energies['A']['acetonitrile'],-energies['B']['acetonitrile']))*(1.602176634e-19*6.02214076e23/1000)
        if abs(actual-state['interaction_energy_kJ_mol_by_medium']['acetonitrile'])>1e-7:
            raise ValueError('Local-contact source energy does not reconstruct')
        groups[state['parent_key']].append(state)
    for key,states in groups.items():
        if sorted(round(s['cap_length_A'],2) for s in states)!=[1.07,1.09,1.11]:
            raise ValueError('Local-contact cap coverage differs')
        center=next(s for s in states if round(s['cap_length_A'],2)==1.09)
        values=[s['interaction_energy_kJ_mol_by_medium']['acetonitrile'] for s in states]
        point=rows[key]
        expected_spacing='adjacent' if key.startswith('occ_2_3') else 'one spacer'
        if (len(point['cap_range_kJ_mol'])!=2 or point['spacing']!=expected_spacing
                or abs(point['centroid_distance_A']-center['parent_contact']['centroid_distance_A'])>1e-9
                or abs(point['interaction_kJ_mol']-center['interaction_energy_kJ_mol_by_medium']['acetonitrile'])>1e-9
                or any(abs(a-b)>1e-9 for a,b in zip(point['cap_range_kJ_mol'],(min(values),max(values))))):
            raise ValueError('Local-contact plotted point differs from source')
    values=[rows['occ_2_3__pose0']['interaction_kJ_mol'],rows['occ_2_4__pose1']['interaction_kJ_mol']]
    if abs(values[0]+24.7615245734)>1e-7 or abs(values[1]+26.2195879462)>1e-7:
        raise ValueError('Local-contact close-pose values differ')
    return values,[panel,metadata,figure,plotter]


def build(evidence,work,review,dipea_evidence,dipea_uv_evidence,tempo_evidence):
    if work.exists():raise FileExistsError('Use a new build directory to preserve previous logs')
    d=json.loads(evidence.read_text())
    for p,h in d['source_sha256'].items():
        if sha(ROOT/p)!=h:raise ValueError('Changed evidence source: '+p)
    anchor=json.loads(dipea_evidence.read_text())
    if (anchor['status']!='primary_DIPEA_salt_NMR_estimate_source_checked_conditional_comparison'
            or anchor['source']['reported_pKa']!=18.7 or anchor['source']['counterion']!='hexafluorophosphate'
            or anchor['calibrated_full_route_prediction']):raise ValueError('Unqualified DIPEA source')
    for path,expected in anchor['source_sha256'].items():
        if sha(path)!=expected:raise ValueError('Changed DIPEA evidence source: '+path)
    uv=json.loads(dipea_uv_evidence.read_text())
    if (uv['status']!='primary_DIPEA_UV_protonation_constant_source_checked_conditional_comparison'
            or uv['source']['reported_value']!=18.64 or uv['source']['temperature_C']!=25.0
            or uv['calibrated_full_route_prediction']):raise ValueError('Unqualified DIPEA UV source')
    for path,expected in uv['source_sha256'].items():
        if sha(path)!=expected:raise ValueError('Changed DIPEA UV evidence source: '+path)
    tempo,tempo_figure=validate_tempo(tempo_evidence)
    contact_values,contact_inputs=validate_local_contacts()
    m=d['benchmark']['metrics']
    values=ROOT/'slides/generated/colleague_values.tex'
    values.write_text('%% Generated from the checked discussion evidence.\n'+
        '\\newcommand{\\BenchmarkMAE}{%.3f}\n\\newcommand{\\BenchmarkRMSE}{%.3f}\n\\newcommand{\\DIPEApKa}{%.1f}\n\\newcommand{\\DIPEAUV}{%.2f}\n\\newcommand{\\LocalAdjacent}{%.2f}\n\\newcommand{\\LocalSpacer}{%.2f}\n'%(m['MAE_V'],m['RMSE_V'],anchor['source']['reported_pKa'],uv['source']['reported_value'],*contact_values))
    inputs=[evidence,dipea_evidence,dipea_uv_evidence,tempo_evidence,tempo_figure,values,Path(__file__).resolve(),ROOT/'slides/redox_workflow.tex',
        ROOT/'slides/generated/colleague_story.tex',ROOT/'slides/generated/colleague_backup.tex',
        ROOT/'results/figures/presentation/figures.json',
        ROOT/'results/figures/polymer_models/pmdi_neighbor_geometry.pdf',*contact_inputs]
    epoch=ROOT/'manuscript/evidence/2026-10-26_native_replay'
    progress_paths=[epoch/n for n in ['polymer_validation.json','finite_ensemble.json','review_package.json']]
    progress=[json.loads(p.read_text()) for p in progress_paths]
    assert progress[0]['qualified_native_G']==progress[0]['expected_graphs']==6
    assert progress[1]['checked_native_G']==progress[1]['expected_native_G']==12
    assert len(progress[1]['finite_pairs'])==6
    assert progress[2]['raw_states']==243 and progress[2]['CPU_Hessian_states']==242
    inputs.extend(progress_paths)
    figures=json.loads((ROOT/'results/figures/presentation/figures.json').read_text())
    plotter=ROOT/'scripts/plotting/presentation/plot_colleague_workflow.py'
    if sha(plotter)!=figures['producer_sha256']:raise ValueError('Changed figure producer')
    inputs.append(plotter)
    for p,h in figures['outputs'].items():
        if sha(p)!=h:raise ValueError('Changed presentation figure: '+p)
        inputs.append(Path(p).resolve())
    assert figures['evidence_sha256']==sha(evidence)
    frozen={str(p.relative_to(ROOT)):sha(p) for p in inputs}
    work.mkdir(parents=True);env={**os.environ,'SOURCE_DATE_EPOCH':'1793037600','FORCE_SOURCE_DATE':'1'}
    logs=[]
    for n in [1,2]:
        log=work/f'pass{n}.log';logs.append(log)
        with log.open('w') as stream:
            subprocess.run(['pdflatex','-interaction=nonstopmode','-halt-on-error','-file-line-error',
                '-recorder','-output-directory',str(work), 'redox_workflow.tex'],cwd=ROOT/'slides',
                stdout=stream,stderr=subprocess.STDOUT,env=env,check=True,timeout=180)
    log=logs[-1].read_text();overfull=re.findall(r'Overfull[^\n]*',log)
    unresolved=[line for line in log.splitlines() if any(k in line for k in ['undefined references','Undefined control sequence','LaTeX Error','Missing character:'])]
    if overfull or unresolved:raise ValueError(dict(overfull=overfull,unresolved=unresolved))
    pdf=work/'redox_workflow.pdf';reader=PdfReader(str(pdf));texts=[p.extract_text() or '' for p in reader.pages]
    if len(reader.pages)!=24:raise ValueError('Expected 17 focused slides and 7 backup slides')
    if '39 / 39' not in texts[7]:raise ValueError('Benchmark slide missing expected coverage')
    if 'Close motifs can interact' not in texts[14]:raise ValueError('Local-contact slide is missing')
    for p,h in frozen.items():
        if sha(ROOT/p)!=h:raise ValueError('Build input changed: '+p)
    for name in ['redox_workflow.pdf','colleague_workflow.pdf']:
        shutil.copy2(pdf,ROOT/'slides'/name)
    row=dict(status='focused_deck_compiler_checked_visual_review_pending',checked=time.time(),
        pages=len(texts),main_slides=17,backup_slides=7,source_sha256=frozen,
        pdf_sha256=sha(pdf),build_directory=str(work),compiler='pdflatex',passes=2,
        overfull_boxes=overfull,unresolved_errors=unresolved,visual_review_performed=False,
        additional_claim_checks=dict(selected_spacing_G=6,finite_pose_G=12,
            finite_pairs=6,reviewer_package_states=243,reviewer_package_Hessians=242,local_contact_triplets=12),
        quantum_jobs_submitted=0,scientific_evidence_modified=False)
    review.mkdir(parents=True,exist_ok=True)
    (review/'build_check.json').write_text(json.dumps(row,indent=2)+'\n')
    (review/'slide_text.json').write_text(json.dumps(texts,indent=2)+'\n')
    print(json.dumps({k:v for k,v in row.items() if k!='source_sha256'},indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--evidence',type=Path,required=True)
    p.add_argument('--dipea-evidence',type=Path,required=True)
    p.add_argument('--dipea-uv-evidence',type=Path,required=True)
    p.add_argument('--tempo-evidence',type=Path,required=True)
    p.add_argument('--work',type=Path,required=True);p.add_argument('--review',type=Path,required=True);a=p.parse_args()
    build(a.evidence.resolve(),a.work.resolve(),a.review.resolve(),a.dipea_evidence.resolve(),a.dipea_uv_evidence.resolve(),a.tempo_evidence.resolve())
