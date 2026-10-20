#!/usr/bin/env python
"""Evidence-based route triage before a new candidate's expensive redox campaign.

Graph validity and a plausible reaction class permit feasibility modeling; neither
is a claim of easy synthesis. Full-property admission needs a reviewed route.
"""
import argparse
import json
from pathlib import Path

from rdkit import Chem

from redox.build.grafting import precursor
from redox.core.common import ROOT
from redox.core.io import atomic_write_json


def screen(row, evidence, assessment):
    result = dict(id=row['id'], route_status='hold', proceed_to_full_redox=False,
                  proceed_to_feasibility_model=False, reasons=[], experimental_ease='unestablished')
    try:
        mol, reaction = precursor(row['frag'])
        result.update(reaction_class=reaction, precursor_smiles=Chem.MolToSmiles(mol),
                      proceed_to_feasibility_model=True)
    except (ValueError, RuntimeError) as exc:
        result.update(route_status='exclude_unreviewed_attachment', reasons=[str(exc)])
        return result
    reviewed = evidence.get(row['id'], {})
    result['evidence'] = reviewed
    hold = assessment.get('candidates', {}).get(row['id'], {}).get('status') == 'on_hold'
    if hold:
        result['reasons'].append('Synthesis selection on hold; realistic base-assisted route unresolved')
    required = ('exact_precursor_access_reviewed', 'attachment_selectivity_reviewed',
                'conditions_defined', 'balanced_route_reviewed', 'redox_core_retention_reviewed')
    missing = [key for key in required if reviewed.get(key) is not True]
    sources=reviewed.get('source_urls')
    if not isinstance(sources,list) or not sources or any(not isinstance(url,str) or not url.startswith('https://') for url in sources):
        missing.append('source_urls')
    conditions=reviewed.get('conditions',{})
    if not isinstance(conditions,dict) or any(not conditions.get(key) for key in ('solvent','reactive_form','base','counterions')):
        missing.append('documented_reaction_conditions')
    result['missing_route_evidence'] = missing
    if reviewed.get('strongly_unfavorable_under_actual_conditions') is True:
        result.update(route_status='exclude_route', proceed_to_feasibility_model=False)
        result['reasons'].append('Reviewed actual-condition attachment route remains strongly unfavorable')
    elif not missing and not hold:
        result.update(route_status='proceed', proceed_to_full_redox=True,
                      experimental_ease='reviewed_workable_route; kinetics/yield not inferred from deltaG alone')
    else:
        result['reasons'].append('Route review required before full redox calculations')
    return result


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--candidates', type=Path, default=ROOT/'results/candidates/current.json')
    ap.add_argument('--evidence', type=Path, default=ROOT/'results/grafting/route_evidence.json')
    ap.add_argument('--out', type=Path, required=True)
    args=ap.parse_args()
    payload=json.loads(args.candidates.read_text())
    rows=payload if isinstance(payload,list) else payload['original39']+payload['additional_candidates']
    evidence=json.loads(args.evidence.read_text()).get('candidates',{}) if args.evidence.exists() else {}
    assessment=json.loads((ROOT/'results/grafting/selection_assessment.json').read_text())
    output=[screen(row,evidence,assessment) for row in rows]
    atomic_write_json(args.out,dict(stage='pre_redox_route_gate',
        scope='Prospective gate; existing 75-candidate calculations are retained',
        candidates=output,counts={status:sum(r['route_status']==status for r in output)
                                 for status in sorted({r['route_status'] for r in output})}))
    print(json.dumps(dict(candidates=len(output),proceed=sum(r['proceed_to_full_redox'] for r in output))))


if __name__=='__main__': main()
