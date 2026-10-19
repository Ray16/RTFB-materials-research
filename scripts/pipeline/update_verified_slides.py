#!/usr/bin/env python
"""Publish the existing Beamer deck only from audited completed candidate data."""
import argparse
import fcntl
import hashlib
import json
import math
import os
import re
from pathlib import Path
import shutil
import signal
import socket
import subprocess
import sys
import time

ROOT=Path(__file__).resolve().parents[2]
FEED=ROOT/'results/candidates/current.json'
LABELS=ROOT/'results/today_2026-10-11_verified_collection/candidate_display_names.json'
SLIDES=ROOT/'slides'
sys.path.insert(0,str(ROOT/'scripts/plotting/candidates'))
from verified_overlay import input_inventory_changed
from collection_layout import artifact_path
from redox.validation.provenance import scoped_spin_summary, input_changes_quarantined, calculation_gid


def save(path,value):
    temp=path.with_name(path.name+f'.{os.getpid()}.tmp')
    temp.write_text(json.dumps(value,indent=2)+'\n');temp.replace(path)


def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def comparison_collection_inputs(payload,root,source):
    """Bind the displayed current figures and include to one completed snapshot."""
    if not any(hook in source.read_text() for hook in (r'\input{generated/candidate_comparison_collection.tex}',r'\input{generated/research_findings.tex}')):return {}
    folder=root/'results/candidates/figures';manifest=folder/'metadata/manifest.json'
    metadata=json.loads(manifest.read_text())
    if metadata['snapshot']!=payload['snapshot_path']:
        raise RuntimeError('Comparison figures require an update from the displayed candidate snapshot')
    inputs={str(manifest):digest(manifest)}
    for name in ('redox_capacity.png','synthesis_tradeoff.png','new_candidate_shifts.png','candidate_structures.png'):
        path=folder/artifact_path(name);expected=metadata['files'][artifact_path(name)]
        if digest(path)!=expected:raise RuntimeError('Comparison figure differs from its checked collection: '+name)
        inputs[str(path)]=expected
    include=root/'slides/generated/candidate_comparison_collection.tex'
    expected=metadata['files']['metadata/candidate_comparison_collection.tex']
    if digest(include)!=expected:raise RuntimeError('Comparison slide include differs from the checked collection')
    inputs[str(include)]=expected
    return inputs


def tex(value):
    mapping={'&':r'\&','%':r'\%','$':r'\$','#':r'\#','_':r'\_',
             '{':r'\{','}':r'\}','~':r'\textasciitilde{}','^':r'\textasciicircum{}',
             '\\':r'\textbackslash{}'}
    return ''.join(mapping.get(c,c) for c in str(value))


def number(value,digits=3):
    return f'{value:.{digits}f}' if type(value) in (int,float) and math.isfinite(value) else '--'


def candidate_label(row,labels):
    cata={'hyp_cata_ndi_ethyl':'Cata-NDI ethyl','hyp_cata_adi_ethyl':'Cata-ADI ethyl'}
    if row['id'] in cata:return cata[row['id']]
    if row['id'].startswith('study_aryl_') and row.get('name'):
        return row['name']
    if labels.get(row['id'],{}).get('short_name'):
        return labels[row['id']]['short_name']
    core=row.get('search_core')
    if core:
        family='Viologen' if core=='viologen' else core.upper()
        contrast=row.get('search_contrast') or ''
        detail=(row.get('search_cap','cap').replace('npropyl','propyl')+' cap'
                if contrast=='cap_length' else contrast.replace('ring_','ring-').replace('linker_C2','ether linker'))
        return family+' '+detail
    for prefix,family in [('hyp_ndi_','NDI'),('hyp_pmdi_','PMDI'),('hyp_ptz_37','PTZ 3,7-')]:
        if row['id'].startswith(prefix):return family+('' if family.endswith('-') else ' ')+row['id'][len(prefix):].replace('2fluoro','2-fluoro')
    return row['id']


def family_screening_frame(rows,when):
    selected=[r for r in rows if r.get('search_core')]
    if not selected:return []
    contrasts={'cap_length':'Propyl cap','ring_F':'Ring F','ring_CN':'Ring CN','ring_OMe':'Ring OMe','ring_Me':'Ring Me','linker_C2':'Ether linker'}
    pages=[]
    for offset in range(0,len(selected),12):
        chunk=selected[offset:offset+12]
        pieces=[r'\begin{frame}{Systematic family screening: comparison panels}',
                r'\scriptsize\centering\renewcommand{\arraystretch}{1.1}',
                r'\begin{tabular}{@{}llrr@{}}\toprule',
                r'Family & Comparison & UMA states & Native states \\ \midrule']
        for r in chunk:
            core='Viologen' if r['search_core']=='viologen' else r['search_core'].upper()
            contrast=contrasts.get(r.get('search_contrast'),r.get('search_contrast') or 'pending')
            pieces.append(tex(core)+' & '+tex(contrast)+' & '+str(r.get('uma_gas_verified_states',0))+'/'+str(r['declared_states'])+
                          ' & '+str(r.get('native_solution_ready_states',0))+'/'+str(r['declared_states'])+r' \\')
        pieces += [r'\bottomrule\end{tabular}\par\vspace{1ex}\raggedright',
                   r'\footnotesize Matched ethyl controls; all ring variants use the same selected site.',
                   r'\par\vspace{1ex}UMA is gas-phase screening. Solution potentials and usable electrons require DFT, thermal, spin and integrity QC.',
                   r'\par\vspace{1ex}Exact synthesis routes and experimental graft performance remain unverified.',
                   f'\\vfill\\scriptsize Audited snapshot: {when}.',r'\end{frame}']
        pages.extend(pieces)
    return pages


def paired_comparison_frame(comparisons,when):
    visible=[r for r in comparisons if r.get('core') in ('viologen','pmdi','ndi') and any(type(r.get(k)) in (int,float) and math.isfinite(r[k])
             for k in ('delta_E1_V','delta_E2_V'))]
    if not visible:return []
    contrasts={'cap_length':'Propyl cap','ring_F':'Ring F','ring_CN':'Ring CN','ring_OMe':'Ring OMe','ring_Me':'Ring Me','linker_C2':'Ether linker'}
    pages=[]
    for offset in range(0,len(visible),10):
        chunk=visible[offset:offset+10]
        pieces=[r'\begin{frame}{Matched family comparisons: admitted potential shifts}',
                r'\scriptsize\centering\renewcommand{\arraystretch}{1.15}',
                r'\begin{tabular}{@{}llrrr@{}}\toprule',
                r'Family & Contrast & $\delta E_1$ & $\delta E_2$ & $\delta(\Delta E)$ \\',
                r'& & \multicolumn{3}{c}{Candidate minus declared control (V)} \\ \midrule']
        for row in chunk:
            core='Viologen' if row['core']=='viologen' else (row['core'] or '').upper()
            pieces.append(tex(core)+' & '+tex(contrasts.get(row['contrast'],row['contrast']))+' & '+
                ' & '.join(number(row.get(k)) for k in ('delta_E1_V','delta_E2_V','delta_wave_spacing_V'))+r' \\')
        pieces += [r'\bottomrule\end{tabular}\par\vspace{1ex}\raggedright',
                   r'\footnotesize Each shift requires an admitted candidate and control event under matched native thermochemistry. Pending events stay blank.',
                   r'\par\vspace{1ex}Positive potential shifts are more positive reductions; the shared Fc reference cancels. A shift alone does not establish improved battery performance.',
                   f'\\vfill\\scriptsize Audited snapshot: {when}.',r'\end{frame}']
        pages.extend(pieces)
    return pages


def expansion_leads_frame(rows,when):
    """Show a diversified follow-up set only when its current ladders are admitted."""
    choices=[('hyp_ptz_37difluoro','3,7-difluoro PTZ'),
             ('search_pmdi_87840c94153a','Ether-linked PMDI'),
             ('search_viologen_1ad45438b9e4','Ring-OMe viologen'),
             ('hyp_pmdi_cyclopropyl','Cyclopropyl PMDI')]
    current={row['id']:row for row in rows}
    selected=[(current[gid],label) for gid,label in choices if gid in current]
    if not selected:return []
    pieces=[r'\begin{frame}{New computational leads for experimental follow-up}',
            r'\small\centering\renewcommand{\arraystretch}{1.25}',
            r'\begin{tabular}{@{}lrrr@{}}\toprule',
            r'Exact graft & $E_1$ & $E_2$ & Nominal capacity \\',
            r'& \multicolumn{2}{c}{V vs Fc/Fc$^+$} & mAh/g \\ \midrule']
    for row,label in selected:
        pieces.append(tex(label)+' & '+number(row.get('E1_vs_Fc_V'))+' & '+
                      number(row.get('E2_vs_Fc_V'))+' & '+
                      number(row.get('capacity_nominal_mAh_g'),1)+r' \\')
    pieces += [r'\bottomrule\end{tabular}\par\vspace{1ex}\raggedright',
               r'\footnotesize Fluoro PTZ is a positive-side voltage hypothesis; its second wave and usable capacity remain provisional under the parent model limitation.',
               r'\par\vspace{1ex}OMe viologen supplies a matched donor-ring contrast; compact PMDI retains useful nominal capacity. Small shifts are unresolved improvements.',
               r'\par\vspace{1ex}This is a follow-up set, not a validated performance ranking. Exact-graft solubility, cycling durability, synthesis and solution kinetics remain unvalidated.',
               f'\\vfill\\scriptsize Corrected native snapshot: {when}.',r'\end{frame}']
    return pieces


def core_topology_rows(payload,identifiers=('design_ndi_ethyl','hyp_cata_ndi_ethyl')):
    rows=payload.get('original39',[])+payload.get('additional_candidates',[])+payload.get('comparison_controls',[])
    by_id={r.get('registry_id') or r['id']:r for r in rows}
    selected=[by_id.get(gid) for gid in identifiers]
    if any(r is None or r.get('solution_potentials_status')!='verified native calculation'
           or r.get('model_primary_admissible') is False
           or any(type(r.get(k)) not in (int,float) or not math.isfinite(r[k])
                  for k in ('E1_vs_Fc_V','E2_vs_Fc_V')) for r in selected):return []
    return selected


def core_topology_frame(payload,when):
    selected=core_topology_rows(payload)
    if not selected:return []
    control,candidate=selected
    pieces=[r'\begin{frame}{Core-topology contrast: cata-NDI and conventional NDI}',
            r'\small\centering\renewcommand{\arraystretch}{1.25}',
            r'\begin{tabular}{@{}lrrr@{}}\toprule',
            r'Exact graft & $E_1$ & $E_2$ & Nominal capacity \\',
            r'& \multicolumn{2}{c}{V vs Fc/Fc$^+$} & mAh/g \\ \midrule']
    for row,label in zip(selected,('Conventional NDI / Et','Cata-NDI / Et')):
        pieces.append(label+' & '+number(row['E1_vs_Fc_V'])+' & '+number(row['E2_vs_Fc_V'])+
                      ' & '+number(row.get('capacity_nominal_mAh_g'),1)+r' \\')
    pieces += [r'\bottomrule\end{tabular}\par\vspace{1.5ex}\raggedright',
        r'\footnotesize Same N-ethyl cap and Merrifield N-graft anchor. This isolates the core-topology isomer contrast.',
        r'\par\vspace{1ex}Candidate minus control: $\delta E_1='+number(candidate['E1_vs_Fc_V']-control['E1_vs_Fc_V'])+
        r'$ V; $\delta E_2='+number(candidate['E2_vs_Fc_V']-control['E2_vs_Fc_V'])+r'$ V. The shared Fc reference cancels.',
        r'\par\vspace{1ex}A more negative anolyte voltage is a lead for follow-up. Overall battery performance, exact synthesis, solubility and solution kinetics remain unvalidated.',
        r'\par\vspace{1ex}\scriptsize Current calculated graft properties; CH$_2$Cl$_2$ source-parent measurements are separate.',
        f'\\vfill\\scriptsize Audited snapshot: {when}.',r'\end{frame}']
    return pieces


def core_extension_frame(payload,when):
    selected=core_topology_rows(payload,('hyp_cata_ndi_ethyl','hyp_cata_adi_ethyl'))
    if not selected:return []
    pieces=[r'\begin{frame}{Fused-ring extension: cata-NDI to cata-ADI}',
            r'\small\centering\renewcommand{\arraystretch}{1.25}',
            r'\begin{tabular}{@{}lrrr@{}}\toprule',
            r'Exact graft & $E_1$ & $E_2$ & Nominal capacity \\',
            r'& \multicolumn{2}{c}{V vs Fc/Fc$^+$} & mAh/g \\ \midrule']
    for row,label in zip(selected,('Cata-NDI / Et','Cata-ADI / Et')):
        pieces.append(label+' & '+number(row['E1_vs_Fc_V'])+' & '+number(row['E2_vs_Fc_V'])+
                      ' & '+number(row.get('capacity_nominal_mAh_g'),1)+r' \\')
    pieces += [r'\bottomrule\end{tabular}\par\vspace{1.5ex}\raggedright',
        r'\footnotesize Same N-ethyl cap and Merrifield anchor. The added fused ring changes both core size and mass.',
        r'\par\vspace{1ex}The voltage/capacity tradeoff is calculated at one solution protocol with field-specific geometry, derivative and spin checks.',
        r'\par\vspace{1ex}The cata-ADI neutral required a reviewed finer CDS derivative step. Its original failed Hessian is preserved; acceptance limits and electronic methods are unchanged.',
        r'\par\vspace{1ex}\scriptsize Exact-graft synthesis, solubility, durability and solution rates remain unvalidated.',
        f'\\vfill\\scriptsize Audited snapshot: {when}.',r'\end{frame}']
    return pieces


def frames(payload,overlays=False):
    summary=payload['summary'];labels=json.loads(LABELS.read_text()) if LABELS.exists() else {}
    rows=[r for r in payload.get('original39',[])+payload.get('additional_candidates',payload.get('extra_NDI',[]))
          if r.get('solution_potentials_status')=='verified native calculation'
          and (r.get('registry_id') or r['id']) not in summary.get('quarantined_input_gids',[])
          and r.get('spin_confidence') not in ('DFT_check_required','wrong_ground_state','pending')
          and type(r.get('E1_vs_Fc_V')) in (int,float) and math.isfinite(r['E1_vs_Fc_V'])]
    if not rows:return None,[]
    rows.sort(key=lambda r:(r.get('family',''),r['E1_vs_Fc_V'],r['id']))
    when=tex(payload['as_of_Chicago'][:19].replace('T',' '))
    pieces=[r'\begin{frame}{Candidate calculations: current verified snapshot}',
        r'\begin{columns}[T]\begin{column}{0.53\textwidth}',
        r'\textbf{39 original proposals + '+str(len(payload.get('additional_candidates',payload.get('extra_NDI',[]))))+r' additional candidates.}',
        r'\begin{itemize}',f"\\item {len(rows)} structures with completed, checked native potential ladders.",
        f"\\item {summary.get('original39_DFT_geometry_states',0)} original-candidate DFT geometry states completed.",
        f"\\item Independent arithmetic: {summary.get('arithmetic_checks',0)} checks, zero discrepancies.",
        r'\item Pending or flagged fields remain blank.',r'\end{itemize}',
        r'\end{column}\begin{column}{0.45\textwidth}',
        r'\begin{block}{What the numbers establish}',
        r'Current geometry, spin evidence, active SMD energies and derivative-checked solution thermochemistry.',
        r'\end{block}',
        r'\footnotesize Predicted molecular properties; experimental accuracy for the exact grafts remains unvalidated.',
        r'\par\vspace{1ex}\scriptsize Native frequency gate: $\nu_{\min}\geq-5$ cm$^{-1}$. Small negative modes within this tolerance remain flagged in the audit.',
        r'\end{column}\end{columns}',
        f'\\vfill\\scriptsize Data snapshot: {when} (America/Chicago).',r'\end{frame}']
    hypotheses=[row for row in payload.get('additional_candidates',[]) if row['id'].startswith('hyp_')]
    if hypotheses:
        pieces += [r'\begin{frame}{Candidate expansion: matched computational hypotheses}',
            r'\scriptsize\centering\renewcommand{\arraystretch}{1.1}',
            r'\begin{tabular}{@{}p{4.8cm}rrr@{}}\toprule',
            r'Design & Nominal capacity & UMA states & Native states \\',
            r'& mAh/g & completed & admitted \\ \midrule']
        for row in hypotheses:
            label=row['id'].replace('hyp_pmdi_','PMDI / ').replace('hyp_ndi_','NDI / ').replace('hyp_ptz_37','PTZ 3,7-')
            if row['id'].startswith('hyp_cata_'):label=candidate_label(row,labels)
            pieces.append(tex(label)+' & '+number(row.get('capacity_nominal_mAh_g'),1)+' & '+
                str(row.get('uma_gas_verified_states',0))+'/'+str(row['declared_states'])+' & '+
                str(row.get('native_solution_ready_states',0))+'/'+str(row['declared_states'])+r' \\')
        pieces += [r'\bottomrule\end{tabular}\par\vspace{1.5ex}\raggedright',
            r'\footnotesize Nominal capacity assumes the declared two-electron ladder. Usable electrons require both potential-window and integrity gates.',
            r'\par\vspace{1ex}UMA is a gas-phase preoptimizer. New solution potentials, disproportionation and solution kinetics remain pending until their separate checks pass.',
            r'\par\vspace{1ex}Exact-graft measurements and synthesis routes are unverified. Parent results are not assigned to these structures.',
            f'\\vfill\\scriptsize Audited snapshot: {when}.',r'\end{frame}']
    pieces += expansion_leads_frame(rows,when)
    pieces += core_topology_frame(payload,when)
    pieces += core_extension_frame(payload,when)
    pieces += family_screening_frame(payload.get('additional_candidates',[]),when)
    pieces += paired_comparison_frame(payload.get('matched_comparisons',[]),when)
    for offset in range(0,len(rows),8):
        chunk=rows[offset:offset+8]
        pieces += [r'\begin{frame}{New candidate potentials with the matched native Fc reference}',
            r'\scriptsize\renewcommand{\arraystretch}{1.25}\centering',
            r'\begin{tabular}{@{}p{4.7cm}rrrrr@{}}\toprule',
            r'Candidate & $E_1$ & $E_2$ & $\Delta E$ & Nominal & Accessible \\',
            r'& \multicolumn{3}{c}{V vs. Fc/Fc$^+$} & \multicolumn{2}{c}{mAh/g} \\ \midrule']
        for r in chunk:
            label=candidate_label(r,labels)
            pieces.append(tex(label)+' & '+' & '.join(number(r.get(k),d) for k,d in (
                ('E1_vs_Fc_V',3),('E2_vs_Fc_V',3),('wave_spacing_V',3),
                ('capacity_nominal_mAh_g',1),('accessible_capacity_mAh_g',1)))+r' \\')
        pieces += [r'\bottomrule\end{tabular}\par',r'\vspace{1ex}\raggedright',
               r'\footnotesize Active-species mass; no counterion/backbone correction. Accessible capacity uses the model window and integrity path; scoped model limits leave it blank.',
            r'\vspace{1ex}\scriptsize Parent residuals do not validate exact grafts. Solution kinetics and synthetic feasibility require separate checks.',
            f'\\vfill\\scriptsize Audited snapshot: {when}.',r'\end{frame}']
    pieces += [r'\begin{frame}{Reference and thermal validation completed}',
        r'\begin{columns}[T]\begin{column}{0.48\textwidth}',
        r'\textbf{Fc geometries repaired on the solution DFT surface.}',
        r'\begin{table}\begin{tabular}{lrr}\toprule',r'State & Earlier & Repaired \\ \midrule',
        'Fc & $-35.1$ & $'+number(summary.get('reference_native_Fc_minimum_frequencies_cm',{}).get('neu'),1)+r'$ \\',
        'Fc$^+$ & $-21.7$ & $'+number(summary.get('reference_native_Fc_minimum_frequencies_cm',{}).get('ox'),1)+r'$ \\',
        r'\bottomrule\end{tabular}\end{table}',
        r'\scriptsize Lowest DFT vibrational frequency (cm$^{-1}$). Both current references are stationary minima.',
        r'\end{column}\begin{column}{0.48\textwidth}\begin{itemize}',
        r'\item Full grid response and matching CPHF grids verified by independent finite gradients.',
        r'\item Matched solution thermochemistry at 298.15 K, including the explicit spin partition.',
        r'\item Five Fe occupation guesses reproduce the chosen geometry-level electronic solution.',
        r'\end{itemize}\end{column}\end{columns}',
        r'\vfill\scriptsize Sources: native thermal certificate, matched Fc reference and occupation-control records (2026-10-11).',
        r'\end{frame}']
    if overlays:
        for title,filename,caption in [
            ('Previous plots with the new candidates: redox and stability','leads_vs_existing_E1_disp_slide.png',
             'Stars and diamonds retain historical table values; native additions are circles. Reference and thermal definitions differ.'),
            ('Previous plots with the new candidates: capacity and gas relaxation','leads_vs_existing_lambda_cap_slide.png',
             'Each new two-electron gas point covers only one event. Gas diagnostics do not validate solution kinetics; missing values are omitted.'),
            ('Previous design space with the new grafted candidates','discovery_capacity_vs_sa_slide.png',
             'New grafts are orange diamonds. Nominal active-species capacity and SA are graph descriptors, with separate QM and route validation.')]:
            image=ROOT/'results/figures/candidates'/filename
            pieces += [r'\begin{frame}{'+title+'}',r'\centering',
                       r'\includegraphics[height=0.67\textheight,width=\linewidth,keepaspectratio]{'+str(image)+'}',
                       r'\vfill\scriptsize '+tex(caption),r'\end{frame}']
    assessment=ROOT/'results/validation/resume_campaign_2026-10-12/candidate_assessment/assessment.json'
    if assessment.exists():
        review=json.loads(assessment.read_text())
        shortlisted=set(review['findings']['practical_shortlist'])
        current={r['id']:r for r in rows}
        pieces += [r'\begin{frame}{Promising computational leads: a diversified test shortlist}',
                   r'\small\renewcommand{\arraystretch}{1.25}',
                   r'\centering\begin{tabular}{@{}lrrr@{}}\toprule',
                   r'Candidate & $E_1$ & $E_2$ & Capacity \\',
                   r'& \multicolumn{2}{c}{V vs Fc/Fc$^+$} & mAh/g \\ \midrule']
        for gid in review['findings']['practical_shortlist']:
            if gid not in current:continue
            r=current[gid]
            label=candidate_label(r,labels)
            pieces.append(tex(label)+' & '+number(r.get('E1_vs_Fc_V'))+' & '+number(r.get('E2_vs_Fc_V'))+' & '+number(r.get('accessible_capacity_mAh_g'),1)+r' \\')
        pieces += [r'\bottomrule\end{tabular}\par',r'\vspace{1ex}\raggedright',
                   r'\footnotesize Compact AQ/diimide/viologen branches plus a conditional high-capacity BQ branch. Capacity is active-species bookkeeping.',
                   r'\vfill\scriptsize Exact-graft experimental performance remains unvalidated; synthesis, solution kinetics and durability are required next gates.',
                   r'\end{frame}',
                   r'\begin{frame}{The trends are plausible; the apparent gains need qualification}',
                   r'\begin{itemize}',
                   r'\item Longer caps add mass and reduce gravimetric capacity; transport benefits have not been demonstrated.',
                   r'\item NDI and PMDI give distinct core/voltage tradeoffs; millivolt alkyl-cap differences are not resolved improvements.',
                   r'\item $\Delta G_{disp}=F\Delta E$: wave spacing and disproportionation reuse the same free-energy information.',
                   r'\item Two potentials do not establish two usable electrons: D3TaLES 05CZHI passes only one step under bond-integrity QC.',
                   r'\end{itemize}',
                   r'\vfill\scriptsize Historical points use superseded reference/thermal definitions. Overall superiority over earlier candidates is unestablished.',
                   r'\end{frame}']
    return '\n'.join(pieces)+'\n',[r['id'] for r in rows]


def publish(payload,out,overlays=False):
    summary=payload['summary']
    if summary.get('arithmetic_discrepancies')!=0 or not input_changes_quarantined(summary):
        return dict(status='waiting_for_audited_snapshot')
    content,ids=frames(payload,overlays)
    if content is None:return dict(status='waiting_for_completed_candidate_potentials')
    # Other candidates continue computing while these completed rows are published.
    # Validate the displayed candidates and their reference, without waiting for the
    # entire fleet's unrelated files to stop changing.
    displayed={r.get('registry_id') or r['id'] for r in
               payload.get('original39',[])+payload.get('additional_candidates',payload.get('extra_NDI',[])) if r['id'] in ids}
    displayed.add('ferrocene')
    registry={r['id']:r.get('registry_id') or r['id'] for r in
              payload.get('original39',[])+payload.get('additional_candidates',[])+payload.get('comparison_controls',[])}
    displayed.update(r.get('registry_id') or r['id'] for r in core_topology_rows(payload))
    displayed.update(r.get('registry_id') or r['id'] for r in core_topology_rows(payload,('hyp_cata_ndi_ethyl','hyp_cata_adi_ethyl')))
    # Paired-shift frames can include a valid first wave from a partial ladder.
    for comparison in payload.get('matched_comparisons',[]):
        if comparison.get('core') in ('viologen','pmdi','ndi') and any(
                type(comparison.get(k)) in (int,float) and math.isfinite(comparison[k])
                for k in ('delta_E1_V','delta_E2_V')):
            displayed.update((registry.get(comparison['id'],comparison['id']),
                              registry.get(comparison['baseline'],comparison['baseline'])))
    summary=scoped_spin_summary(summary,displayed,ROOT)
    if input_inventory_changed(summary):return dict(status='waiting_for_current_input_inventory')
    inputs=dict(summary.get('source_sha256',{}))
    for p,h in summary.get('admitted_input_sha256',{}).items():
        owner=calculation_gid(p,ROOT)
        if owner is None or owner in displayed:inputs[p]=h
    if any(not Path(p).exists() or digest(p)!=h for p,h in inputs.items()):
        return dict(status='waiting_for_current_snapshot_inputs')
    source=SLIDES/'redox_workflow.tex';source_hash=digest(source)
    if overlays:
        overlay_feed=out/f'overlay_input_{time.time_ns()}.json'
        save(overlay_feed,payload)
        argv=[sys.executable,str(ROOT/'scripts/plotting/candidates/update_candidate_overlays.py'),
              '--feed',str(overlay_feed),'--status',str(out/'overlay_status.json')]
        with (out/'overlay.log').open('w') as output:
            proc=subprocess.run(argv,cwd=ROOT,stdout=output,stderr=subprocess.STDOUT,timeout=180)
        if proc.returncode:raise RuntimeError('Existing plot overlay update failed; previous PDF preserved')
    inputs.update(comparison_collection_inputs(payload,ROOT,source))
    if r'\InputIfFileExists{generated/verified_candidates.tex}' not in source.read_text():
        raise RuntimeError('Verified-slide hook is missing; preserve user edits and inspect')
    generated=SLIDES/'generated';generated.mkdir(exist_ok=True)
    include=generated/'verified_candidates.tex';temp=include.with_name('.verified_candidates.tmp')
    temp.write_text(content);temp.replace(include)
    main_frames=None
    if r'\input{generated/research_findings.tex}' in source.read_text():
        from research_presentation import research_frames
        main_content=research_frames(payload,ROOT)
        main_include=generated/'research_findings.tex'
        main_include.write_text(main_content)
        main_frames=len(re.findall(r'\\begin\{frame\}',main_content))
        for path in (ROOT/'scripts/plotting/candidates/research_presentation.py',
                     ROOT/'results/candidates/figures/data/highlights.json',
                     ROOT/'results/discovery/vinyl_monomers/proposals.json'):
            inputs[str(path)]=digest(path)
    # Bind supporting TeX as well as the main story to the checked build.
    for path in list(generated.glob('*.tex'))+list((SLIDES/'backup').glob('*.tex')):
        inputs[str(path)]=digest(path)
    build=out/'builds'/str(time.time_ns());build.mkdir(parents=True)
    argv=['pdflatex','-interaction=nonstopmode','-halt-on-error',f'-output-directory={build}',str(source)]
    for pass_number in (1,2):
        with (build/f'pass{pass_number}.log').open('w') as log:
            proc=subprocess.run(argv,cwd=SLIDES,stdout=log,stderr=subprocess.STDOUT,timeout=120)
        if proc.returncode:raise RuntimeError(f'Slide build failed; prior published PDF preserved: {build}')
    pdf=build/'redox_workflow.pdf'
    from pypdf import PdfReader
    document=PdfReader(pdf)
    text=' '.join(' '.join((page.extract_text() or '').split()) for page in document.pages)
    pages=len(document.pages)
    if pages<16 or 'New candidate potentials' not in text:
        raise RuntimeError('Built PDF is missing the generated result frames')
    if main_frames is not None:
        if document.pages[main_frames].extract_text().strip().splitlines()[0]!='Backup':
            raise RuntimeError('Supporting slides must follow the main research story in Backup')
    if digest(source)!=source_hash or input_inventory_changed(summary) or any(not Path(p).exists() or digest(p)!=h for p,h in inputs.items()):
        return dict(status='source_changed_during_build',preserved_build=str(build))
    shutil.copy2(pdf,SLIDES/'.redox_workflow.checked.pdf')
    (SLIDES/'.redox_workflow.checked.pdf').replace(SLIDES/'redox_workflow.pdf')
    return dict(status='published_checked_pdf',snapshot=payload['snapshot_path'],candidate_ids=ids,
                slides=str(SLIDES/'redox_workflow.pdf'),main_research_slides=main_frames,backup_slides=(pages-main_frames if main_frames is not None else None),pages=pages,pdf_sha256=digest(SLIDES/'redox_workflow.pdf'),
                source_sha256=source_hash,build=str(build),updated=time.time())


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--out',type=Path,required=True)
    ap.add_argument('--watch',action='store_true');ap.add_argument('--interval',type=int,default=120)
    ap.add_argument('--feed',type=Path,default=FEED)
    ap.add_argument('--max-hours',type=float,help='optional watch duration spanning midnight')
    ap.add_argument('--update-plots',action='store_true',help='add expansion points to existing comparison plots')
    a=ap.parse_args();out=a.out.resolve();out.mkdir(parents=True,exist_ok=True)
    lock=(out/'publisher.lock').open('a+');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    status=dict(host=socket.gethostname(),pid=os.getpid(),cpu_affinity=sorted(os.sched_getaffinity(0)),started=time.time())
    from redox.workflows.publication import watch_deadline
    try:deadline=watch_deadline(a.max_hours)
    except ValueError as exc:ap.error(str(exc))
    last=None
    def interrupted(signum,frame):raise KeyboardInterrupt(f'signal {signum}')
    signal.signal(signal.SIGTERM,interrupted)
    try:
        while True:
            if a.feed.exists():
                payload=json.loads(a.feed.read_text())
                if payload['snapshot_path']!=last:
                    result=publish(payload,out,a.update_plots);status.update(result);save(out/'status.json',status)
                    if result['status']=='published_checked_pdf':last=payload['snapshot_path']
            if not a.watch or time.time()>=deadline:break
            time.sleep(a.interval)
    except (Exception,KeyboardInterrupt) as exc:
        status.update(status='stopped_preserving_published_pdf',error=repr(exc));save(out/'status.json',status);raise
    finally:lock.close()


if __name__=='__main__':main()
