#!/usr/bin/env python
"""Assess a reported DIPEAH/PF6 experimental pKa against frozen model cycles."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import re
import xml.etree.ElementTree as ET

ROOT=Path(__file__).resolve().parents[2]
EPOCH=ROOT/'manuscript/evidence/2026-10-26_native_replay'
R=8.31446261815324


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def assess(article,base_path,acid_path,origin_path):
    root=ET.parse(article).getroot()
    dois=[node.text.lower() for node in root.findall('.//article-meta/article-id') if node.get('pub-id-type')=='doi']
    if dois!=['10.1039/d3sc02611f']:raise ValueError('Wrong primary article identity')
    title=' '.join(''.join(root.find('.//article-meta/title-group/article-title').itertext()).split())
    if not title.startswith('Controlling product selectivity during dioxygen reduction'):
        raise ValueError('Wrong primary article title')
    paragraphs=[' '.join(' '.join(p.itertext()).split()) for p in root.iter('p')]
    relevant=[p for p in paragraphs if '18.7' in p and 'DIPEAHPF' in p and 'competition experiments' in p]
    if len(relevant)!=1 or 'MeCN' not in relevant[0] or 'hexafluorophosphate' not in relevant[0]:
        raise ValueError('Primary NMR measurement, medium and salt identity not established')
    base=json.loads(base_path.read_text());acid=json.loads(acid_path.read_text());origin=json.loads(origin_path.read_text())
    hcl=acid['acid_anchor']['pKa_MeCN'];temperature=acid['temperature_K']
    dbu=next(r['reference']['pKaH_MeCN'] for r in acid['rows'] if r['base']=='dbu')
    if hcl!=10.30 or dbu!=24.31 or temperature!=298.15:raise ValueError('Reference condition changed')
    dft=base['corrections']['dipea']['free_ions_kJ_mol']
    dft_dbu=base['corrections']['dbu']['free_ions_kJ_mol']
    relative_dft=origin['DFT_relative_DIPEA_DBU_transfer_kJ_mol']
    if abs(relative_dft-(dft-dft_dbu))>.001:raise ValueError('Relative base transfer does not match frozen components')
    nmr=18.7;factor=R*temperature*math.log(10)/1000
    reference_capture=-factor*(nmr-hcl)
    reference_relative=factor*(dbu-nmr)
    if not all(math.isfinite(v) for v in (reference_capture,reference_relative,relative_dft,dft)):
        raise ValueError('Nonfinite comparison')
    candidates=[]
    for row in base['candidates']:
        reference=row['reference_deltaG_kJ_mol'];primary=row['base_models']['dipea']['free_ions_kJ_mol']
        if reference is None or primary is None:
            if reference is not None or primary is not None:raise ValueError('Partial DIPEA sum is inconsistent')
            candidates.append(dict(id=row['id'],reference_deltaG_kJ_mol=None,all_DFT_DIPEA_kJ_mol=None,
                hybrid_capture_sensitivity_kJ_mol=None,hybrid_is_calibrated_reaction_prediction=False))
            continue
        if abs(reference+dft-primary)>1e-7:raise ValueError('Primary DIPEA sum changed')
        candidates.append(dict(id=row['id'],reference_deltaG_kJ_mol=reference,
            all_DFT_DIPEA_kJ_mol=primary,hybrid_capture_sensitivity_kJ_mol=reference+reference_capture,
            hybrid_is_calibrated_reaction_prediction=False))
    source_paths=[article,base_path,acid_path,origin_path,Path(__file__).resolve()]
    return dict(status='primary_DIPEA_salt_NMR_estimate_source_checked_conditional_comparison',
        source=dict(doi=dois[0],title=title,article_url='https://doi.org/10.1039/D3SC02611F',
            open_fulltext_url='https://www.ebi.ac.uk/europepmc/webservices/rest/PMC10952101/fullTextXML',
            publication='Chemical Science 2024, 15, 4478–4488',
            authors=['Emma N. Cook','Ian M. Courter','Diane A. Dickie','Charles W. Machan'],
            observation='1H-NMR titration/competition estimate for diisopropylethylammonium hexafluorophosphate',
            solvent='acetonitrile',reported_pKa=18.7,counterion='hexafluorophosphate',
            source_temperature_and_activity_conventions_verified=False,
            locator='Experimental: synthesis of diisopropylethylammonium hexafluorophosphate',
            article_license='CC BY-NC 3.0; raw XML retained on shared storage'),
        comparison=dict(temperature_K=temperature,R_J_mol_K=R,RTln10_kJ_mol=factor,
            HCl_pKa=hcl,DBUH_pKa=dbu,DIPEAH_reported_pKa=nmr,
            DFT_separated_ion_capture_kJ_mol=dft,pKa_implied_capture_kJ_mol=reference_capture,
            DFT_minus_pKa_implied_capture_kJ_mol=dft-reference_capture,
            DFT_relative_DIPEA_DBU_transfer_kJ_mol=relative_dft,
            pKa_implied_DIPEA_DBU_transfer_kJ_mol=reference_relative,
            DFT_minus_pKa_implied_relative_kJ_mol=relative_dft-reference_relative,
            relative_transfer_reaction='DIPEA + DBUH+ -> DIPEAH+ + DBU',
            separation_ion_and_activity_conventions_equivalent=False),
        candidate_hybrid_sensitivities=candidates,
        limitations=['Experimental number is a salt-specific NMR estimate; PF6 association and measurement conditions may affect comparison.',
            'HCl value is an author-tabulated recommendation, not necessarily matched in temperature or activity convention.',
            'HCl cancels in the full DFT grafting cycle; substituting only the capture term changes error cancellation.',
            'Agreement of relative-base transfer does not uniquely identify the cause of the HCl capture discrepancy.',
            'No experimental grafting free energy, selectivity, yield or material admission is established.'],
        calibrated_full_route_prediction=False,primary_grafting_results_modified=False,candidate_decisions_modified=False,
        source_sha256={str(p):sha(p) for p in source_paths})


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--article',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    if a.out.exists():raise FileExistsError('Preserve prior reference assessments')
    result=assess(a.article.resolve(),EPOCH/'base_assisted.json',EPOCH/'acid_base.json',EPOCH/'basicity_reference_origin.json')
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result['comparison'],indent=2))


if __name__=='__main__':main()
