#!/usr/bin/env python
"""Compare a 25 °C DIPEA spectrophotometric protonation measurement to DFT."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import re
import xml.etree.ElementTree as ET

from pypdf import PdfReader

ROOT=Path(__file__).resolve().parents[2]
EPOCH=ROOT/'manuscript/evidence/2026-10-26_dipea_nmr'


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def extract_observation(article_path,si_path):
    article=ET.parse(article_path).getroot()
    dois=[x.text.lower() for x in article.findall('.//article-meta/article-id') if x.get('pub-id-type')=='doi']
    if dois!=['10.1039/d1ra04738h']:raise ValueError('Wrong experimental source DOI')
    paragraphs=[' '.join(' '.join(p.itertext()).split()) for p in article.iter('p')]
    text=next((p for p in paragraphs if 'Protonation constants of DIPEA in MeCN' in p and 'bromocresol green' in p.lower()),None)
    if text is None or '25' not in text:raise ValueError('Missing DIPEA MeCN measurement and indicator details')
    pdf=PdfReader(str(si_path));page=' '.join((pdf.pages[42].extract_text() or '').split())
    if not ('Figure S51.' in page and 'DIPEA' in page and 'MeCN' in page
            and '25.0' in page and re.search(r'log\s+K\s*H\s*=\s*18\.64\(1\)',page)):
        raise ValueError('SI page does not establish the reported MeCN protonation constant')
    return dict(doi=dois[0],article_url='https://doi.org/10.1039/D1RA04738H',
        title='Protonation and anion-binding properties of aromatic sulfonylurea derivatives',
        publication='RSC Advances 2021, 11, 23992–24000',
        measured_quantity='log K^H for DIPEA protonation',reported_value=18.64,reported_notation='18.64(1)',
        solvent='acetonitrile',temperature_C=25.0,temperature_variation_C=0.1,
        technique='spectrophotometric titration with bromocresol green',
        fixed_indicator_logK_values=[18.5,11.0],measurement_locator='Supporting Information Figure S51, printed and PDF page 43',
        comparison_conditions_fully_matched=False)


def assess(article,si,base_path,acid_path,origin_path,nmr_path):
    source=extract_observation(article,si);base=json.loads(base_path.read_text());acid=json.loads(acid_path.read_text())
    origin=json.loads(origin_path.read_text());nmr=json.loads(nmr_path.read_text())
    if nmr['source']['reported_pKa']!=18.7 or nmr['source']['counterion']!='hexafluorophosphate':
        raise ValueError('NMR source identity changed')
    hcl=acid['acid_anchor']['pKa_MeCN'];t=acid['temperature_K'];dbu=next(r['reference']['pKaH_MeCN'] for r in acid['rows'] if r['base']=='dbu')
    if (hcl,t,dbu)!=(10.30,298.15,24.31):raise ValueError('Reference conditions changed')
    dft=base['corrections']['dipea']['free_ions_kJ_mol'];dft_relative=origin['DFT_relative_DIPEA_DBU_transfer_kJ_mol']
    factor=8.31446261815324*t*math.log(10)/1000;uv=source['reported_value']
    capture=-factor*(uv-hcl);relative=factor*(dbu-uv)
    rows=[]
    for r in base['candidates']:
        old=r['reference_deltaG_kJ_mol'];primary=r['base_models']['dipea']['free_ions_kJ_mol']
        if old is None or primary is None:
            if old is not None or primary is not None:raise ValueError('Inconsistent missing primary reaction')
            value=None
        else:
            if abs(old+dft-primary)>1e-7:raise ValueError('Primary DIPEA reaction does not close')
            value=old+capture
        rows.append(dict(id=r['id'],reference_deltaG_kJ_mol=old,all_DFT_DIPEA_kJ_mol=primary,
                         hybrid_capture_sensitivity_kJ_mol=value,calibrated_reaction_prediction=False))
    inputs=[article,si,base_path,acid_path,origin_path,nmr_path,Path(__file__).resolve()]
    return dict(status='primary_DIPEA_UV_protonation_constant_source_checked_conditional_comparison',source=source,
        independent_NMR_salt_estimate_pKa=18.7,comparison=dict(temperature_K=t,HCl_pKa=hcl,DBUH_pKa=dbu,
            RTln10_kJ_mol=factor,DFT_capture_kJ_mol=dft,pKa_implied_capture_kJ_mol=capture,
            DFT_minus_pKa_implied_capture_kJ_mol=dft-capture,DFT_relative_base_transfer_kJ_mol=dft_relative,
            pKa_implied_relative_base_transfer_kJ_mol=relative,DFT_minus_pKa_implied_relative_kJ_mol=dft_relative-relative,
            relative_reaction='DIPEA + DBUH+ -> DIPEAH+ + DBU'),candidate_hybrid_sensitivities=rows,
        limitations=['Measured protonation constant is anchored to a bromocresol-green indicator scale.',
            'Exact activity and ionic-association conventions are not established to match the separate HCl and DBU sources.',
            'The 2024 NMR estimate uses a PF6 salt and MeCN-d3; neither measurement directly observes chloride pairing.',
            'HCl cancels from the full all-DFT reaction; hybrid substitutions change error cancellation and do not calibrate the whole route.'],
        calibrated_full_route_prediction=False,primary_grafting_results_modified=False,candidate_decisions_modified=False,
        source_sha256={str(p):sha(p) for p in inputs})


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--article',type=Path,required=True);p.add_argument('--si',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    if a.out.exists():raise FileExistsError('Preserve existing experimental assessment')
    result=assess(a.article.resolve(),a.si.resolve(),EPOCH/'base_assisted.json',EPOCH/'acid_base.json',EPOCH/'basicity_reference_origin.json',EPOCH/'dipea_experimental_anchor.json')
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result['comparison'],indent=2))
