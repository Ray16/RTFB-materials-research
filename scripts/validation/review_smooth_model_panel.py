#!/usr/bin/env python
"""Join component and full quantum controls for a neutral CHNOCl model pilot."""
import argparse
import hashlib
import json
from pathlib import Path
import time

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api


def checked_controls(root):
    path=Path(root)/'qualification.json';row=json.loads(path.read_text())
    if (row['status']!='experimental_neutral_CHNOCl_model_controls_qualified'
            or row['primary_candidate_admission'] or row['full_molecular_polymer_G_qualified']):
        raise ValueError('Qualified neutral component/quantum model panel required')
    api.verify_sources(row['source_sha256'])
    return {str(path):api.digest(path),**row['source_sha256']}


def review(panels,gradient_panel,out):
    from ase.io import read
    if len(panels)!=2 or len(set(panels))!=2:raise ValueError('Two separate two-spin panels required')
    checker=api.load_module('four_chemical_control_panel',ROOT/'scripts/validation/review_smooth_quantum.py')
    model=api.load_module('smooth_panel_model_identity',ROOT/'scripts/pipeline/smooth_cds_control.py')
    proof={};coverage=set();controls=[]
    identity=hashlib.sha256(json.dumps(model.MODEL,sort_keys=True).encode()).hexdigest()
    for panel in panels:
        proof.update(checker.checked_controls(panel))
        report=json.loads((panel/'independent_review.json').read_text())
        if report['portable_method_id']!=identity:raise ValueError('Different surface models cannot be joined')
        for mult in (1,2):
            path=panel/'controls'/str(mult)/'status.json';row=json.loads(path.read_text())
            xyz=Path(row['raw_root'])/'input.xyz';atoms=read(xyz);coverage.update(atoms.get_chemical_symbols())
            controls.append(dict(panel=str(panel),mult=mult,atoms=atoms.get_chemical_symbols(),status_path=str(path)))
    if coverage!={'C','H','N','O','Cl'}:raise ValueError('Requested neutral element coverage is incomplete')
    gradients=json.loads(gradient_panel.read_text())
    if gradients['status']!='smooth_amine_imide_full_gradients_checked' or set(gradients['states'])!={'precursor','base','protonated','pair'}:
        raise ValueError('Full imide/amine component-gradient checks required')
    for row in gradients['states'].values():
        if row['tested_components']!=row['total_components'] or row['maximum_error_Ha_Bohr']>1e-6:
            raise ValueError('Full Cartesian gradient panel did not pass')
    api.verify_sources(gradients['source_sha256']);proof.update(gradients['source_sha256'])
    proof.update({str(p):api.digest(p) for p in (Path(__file__).resolve(),gradient_panel,ROOT/'scripts/pipeline/smooth_cds_control.py')})
    out.mkdir(parents=True,exist_ok=False)
    result=dict(status='experimental_neutral_CHNOCl_model_controls_qualified',checked=time.time(),
        method=model.MODEL,portable_method_id=identity,quantum_controls=controls,
        full_component_gradients=list(gradients['states']),source_sha256=proof,
        primary_candidate_admission=False,full_molecular_polymer_G_qualified=False,
        scope='Four small neutral full-quantum controls plus component gradients, supporting a separately gated large neutral geometry pilot. Charged molecular QM and polymer thermochemistry remain unqualified.')
    atomic_write_json(out/'qualification.json',result)
    print(json.dumps(dict(status=result['status'],small_quantum_controls=len(controls))))


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--panel',type=Path,action='append',required=True)
    parser.add_argument('--gradient-panel',type=Path,required=True);parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args();review([p.resolve() for p in args.panel],args.gradient_panel.resolve(),args.out.resolve())


if __name__=='__main__':main()
