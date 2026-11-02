#!/usr/bin/env python
"""Independent source and arithmetic review of the DIPEA UV anchor."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import xml.etree.ElementTree as ET

from pypdf import PdfReader


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def review(summary_path,article,si,nmr_si,out):
    s=json.loads(summary_path.read_text())
    for name,expected in s['source_sha256'].items():
        if sha(name)!=expected:raise ValueError('Source changed: '+name)
    doi=[x.text.lower() for x in ET.parse(article).getroot().findall('.//article-meta/article-id') if x.get('pub-id-type')=='doi']
    page=' '.join((PdfReader(str(si)).pages[42].extract_text() or '').split())
    if doi!=['10.1039/d1ra04738h'] or not all(x in page for x in ('Figure S51.','MeCN','25.0','log KH = 18.64(1)')):
        raise ValueError('Independent check of the experimental source failed')
    nmr=' '.join((PdfReader(str(nmr_si)).pages[12].extract_text() or '').split())
    if not all(x in nmr for x in ('18.81','0.02 M','MeCN-d3','0.142','0.083')):
        raise ValueError('NMR supporting-information context not verified')
    c=s['comparison'];f=8.31446261815324*298.15*math.log(10)/1000
    expected_capture=-f*(18.64-10.30);expected_relative=f*(24.31-18.64)
    if (abs(c['pKa_implied_capture_kJ_mol']-expected_capture)>1e-10
            or abs(c['pKa_implied_relative_base_transfer_kJ_mol']-expected_relative)>1e-10
            or abs(c['DFT_minus_pKa_implied_capture_kJ_mol']-(c['DFT_capture_kJ_mol']-expected_capture))>1e-10
            or abs(c['DFT_minus_pKa_implied_relative_kJ_mol']-(c['DFT_relative_base_transfer_kJ_mol']-expected_relative))>1e-10
            or s['calibrated_full_route_prediction'] or s['candidate_decisions_modified']):
        raise ValueError('Experimental comparison does not reconstruct')
    inputs=[summary_path,article,si,nmr_si,Path(__file__).resolve()]
    report=dict(status='independent_UV_source_and_comparison_review_checked',
        source_uv_doi=doi[0],UV_reported_logKH='18.64(1)',UV_temperature_C=25.0,
        UV_indicator='bromocresol green',NMR_reported_pKa=18.7,NMR_medium='MeCN-d3',
        NMR_initial_salt_molar=0.02,NMR_reference_pKa=18.81,
        NMR_SI_reports_absolute_pKa_difference=0.142,NMR_SI_reported_SD=0.083,
        NMR_table_sign_convention='Table S3 reports negative entries under a Ka-ratio heading and Eq. S3 uses an absolute value; use the article-reported 18.7 without an independent signed re-fit.',
        direct_chloride_equivalence_inferred=False,full_grafting_accuracy_inferred=False,
        DFT_minus_UV_implied_capture_kJ_mol=c['DFT_minus_pKa_implied_capture_kJ_mol'],
        DFT_minus_UV_implied_relative_kJ_mol=c['DFT_minus_pKa_implied_relative_kJ_mol'],
        source_sha256={str(p):sha(p) for p in inputs})
    if out.exists():raise FileExistsError('Preserve existing source review')
    out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(report,indent=2)+'\n');return report


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--summary',type=Path,required=True);p.add_argument('--article',type=Path,required=True);p.add_argument('--si',type=Path,required=True);p.add_argument('--nmr-si',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    r=review(a.summary.resolve(),a.article.resolve(),a.si.resolve(),a.nmr_si.resolve(),a.out.resolve());print(json.dumps({k:v for k,v in r.items() if k!='source_sha256'},indent=2))
