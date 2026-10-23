#!/usr/bin/env python
"""Keep computed/reference-origin ambiguities out of experimental validation."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import re
import time
import unicodedata


def origin_role(record, purpose):
    measured=record["origin_class"] in ("experimental_measurement","author_recommended_measured_scale")
    eligible=measured and not record.get("origin_conflict",False)
    if record.get("independent_experimental_anchor") is True and not eligible:
        raise ValueError("Computed or ambiguous source was promoted to an experimental anchor")
    if purpose=="experimental_validation" and not eligible:
        raise ValueError("An independent measured source is required for experimental validation")
    if purpose not in ("experimental_validation","model_sensitivity"):
        raise ValueError("Unknown reference purpose")
    return dict(independent_experimental_origin_eligible=eligible,
        full_route_accuracy_validated=False,role=purpose)


def reviewed_DIPEA_context(record):
    if (record['base']!='dipea' or record['solvent']!='acetonitrile' or record['value_type']!='pKaH'
            or record['reported_value']!=18.6 or record['origin_class']!='computed_PCM_estimate'
            or not record['origin_conflict']):
        raise ValueError('Reference value/context differs from the source-specific DIPEA review')


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args()
    from redox.core.common import ROOT
    from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity
    from redox.workflows.checkpoints import verify_state_checkpoint
    from rdkit import Chem
    from pypdf import PdfReader
    import os
    with allocation_guard():host_capacity(1);os.sched_setaffinity(0,cpu_allocation(1))
    args.out.mkdir(parents=True,exist_ok=False)
    spec_path=ROOT/'config/basicity_reference_origins.json';base_path=ROOT/'results/grafting/base_assisted/current.json'
    acid_path=ROOT/'results/validation/acid_base_external_2026-10-22/comparison.json'
    spec=json.loads(spec_path.read_text())['references'][0];base=json.loads(base_path.read_text());acid=json.loads(acid_path.read_text())
    reviewed_DIPEA_context(spec)
    pdf=ROOT/spec['raw_pdf'];sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
    if sha(pdf)!=spec['raw_pdf_sha256']:raise ValueError('Reviewed source PDF changed')
    pages=PdfReader(pdf).pages
    normalize=lambda s:re.sub('[^a-z0-9]','',unicodedata.normalize('NFKD',s).encode('ascii','ignore').decode().lower())
    if ('estimatedusingthepcmmodel' not in normalize(pages[1].extract_text())
            or 'experimentalpkah' not in normalize(pages[2].extract_text())
            or spec['origin_class']!='computed_PCM_estimate' or not spec['origin_conflict']):
        raise ValueError('Specific-note/general-label review does not match the inspected source')
    if Chem.MolToInchiKey(Chem.MolFromSmiles(spec['smiles']))!=Chem.MolToInchiKey(Chem.MolFromSmiles(base['bases']['dipea']['neutral'])):
        raise ValueError('Literature base identity differs from the computed base')
    role=origin_role(spec,'model_sensitivity')
    components={};sources={str(p):sha(p) for p in (Path(__file__).resolve(),spec_path,base_path,acid_path,pdf)}
    campaign=ROOT/'results/grafting/campaigns/2026-10-19_base_assisted_pilot'
    for label in ('dbu_neutral','dbu_protonated','dipea_neutral','dipea_protonated'):
        path=campaign/'done'/(label+'.json');report=verify_state_checkpoint(path);components[label]=report['G_eV'];sources[str(path)]=sha(path)
    temperature=298.15;rtln10=.00831446261815324*temperature*math.log(10)
    capture=-rtln10*(spec['reported_value']-acid['acid_anchor']['pKa_MeCN'])
    relative=(components['dipea_protonated']+components['dbu_neutral']-components['dipea_neutral']-components['dbu_protonated'])*96.48533212331002
    dbu_reference=next(r['reference']['pKaH_MeCN'] for r in acid['rows'] if r['base']=='dbu')
    inferred=dbu_reference-relative/rtln10
    sensitivity=[]
    for row in base['candidates']:
        if row['reference_deltaG_kJ_mol'] is None:continue
        sensitivity.append(dict(id=row['id'],reference_deltaG_kJ_mol=row['reference_deltaG_kJ_mol'],
            all_DFT_DIPEA_kJ_mol=row['base_models']['dipea']['free_ions_kJ_mol'],
            fixed_reference_plus_computed_literature_capture_kJ_mol=row['reference_deltaG_kJ_mol']+capture,
            calibrated_full_route_prediction=False))
    result=dict(status='source_origin_reviewed_no_independent_DIPEA_experiment_admitted',checked=time.time(),source=spec,
        source_role=role,independent_experimental_pKaH_MeCN=None,comparison_temperature_K=temperature,
        source_temperature_verified=False,HCl_cancels_in_complete_all_DFT_route=True,
        computed_literature_capture_kJ_mol=capture,DFT_relative_DIPEA_DBU_transfer_kJ_mol=relative,
        DIPEA_pKaH_inferred_from_DFT_relative_transfer_and_DBU_anchor=inferred,
        comparison_with_literature_is_model_to_model=True,candidate_model_sensitivities=sensitivity,
        primary_grafting_results_modified=False,candidate_decisions_modified=False,source_sha256=sources)
    (args.out/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
    (args.out/'source_spec.json').write_bytes(spec_path.read_bytes())
    print(json.dumps({k:v for k,v in result.items() if k not in ('source','source_sha256','candidate_model_sensitivities')}))


if __name__=='__main__':main()
