#!/usr/bin/env python
"""Fresh local single points in an independently measured execution class.

The complete original coordinates are copied for each selected input. Only the
execution class transfers; no reference energy, orbitals or density is reused.
"""
import argparse
import json
import os
from pathlib import Path
import time

from rdkit import Chem

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api


def prepare(construction,resource_review,selections,out):
    witness=json.loads(resource_review.read_text());api.verify_sources(witness['source_sha256'])
    if witness['status']!='frozen_local_electronic_state_and_resource_pilot_independently_checked':
        raise ValueError('Completed independent local resource evidence required')
    point=api.read_plan(Path(witness['campaign']))
    worker=api.load_module('local_electronic_unchanged_calculator',ROOT/'scripts/pipeline/local_motif_energy.py')
    if witness['settings']!=worker.SETTINGS or point['settings']!=worker.SETTINGS:
        raise ValueError('Local electronic protocol differs from the measured point')
    manifest=json.loads((construction/'manifest.json').read_text());api.verify_sources(manifest['source_sha256'])
    structural=ROOT/'results/reviews/2026-10-26/local_motif_construction/summary.json'
    checked=json.loads(structural.read_text());api.verify_sources(checked['source_sha256'])
    if checked['source_sha256'].get(str(construction/'manifest.json'))!=api.digest(construction/'manifest.json'):
        raise ValueError('Construction review belongs to another input set')
    out.mkdir(parents=True,exist_ok=False);entries=[]
    expected=Chem.MolToSmiles(Chem.MolFromSmiles('CN1C(=O)c2ccccc2C1=O'))
    for selection in selections:
        key,label=selection.rsplit(':',1);row=next(r for r in manifest['entries'] if r['key']==key)
        item=next(r for r in row['molecules'] if r['id']==label)
        molecule=Chem.MolFromMolFile(item['molecular_graph'],removeHs=False)
        if (item['atoms']!=point['task']['atoms'] or item['charge']!=0 or molecule is None
                or Chem.GetFormalCharge(molecule)!=0
                or any(Chem.MolToSmiles(Chem.RemoveHs(m))!=expected for m in Chem.GetMolFrags(molecule,asMols=True))):
            raise ValueError('Local molecule lies outside the measured constitution/state/size class')
        name=key+'__'+label;case=out/name;case.mkdir();xyz=case/'input.xyz';xyz.write_bytes(Path(item['geometry']).read_bytes())
        inputs={**witness['source_sha256'],**point['input_sha256']}
        for path in (construction/'manifest.json',structural,resource_review,Path(item['geometry']),Path(item['molecular_graph']),xyz):inputs[str(path)]=api.digest(path)
        plan=dict(point,task=dict(point['task'],key=name,source_parent=row['parent_key']),input_xyz=str(xyz),
            geometry_sha256=api.digest(xyz),input_sha256=inputs,
            producer_sha256={**point['producer_sha256'],str(Path(__file__).resolve()):api.digest(__file__)},
            raw_root=str(ROOT/'calcs/local_motif_interactions'/out.name/name),created=time.time(),
            resource_reference=witness['campaign'],resource_review=str(resource_review),source_input_key=key,source_input_label=label,
            reference_energy_transferred=False,prospective_resource_pilot=False)
        atomic_write_json(case/'plan.json',plan);entries.append(dict(key=name,campaign=str(case),plan_sha256=api.digest(case/'plan.json')))
    atomic_write_json(out/'campaign.json',dict(status='same_class_fresh_local_energies_prepared',entries=entries,
        source_construction=str(construction),source_review=str(structural),resource_review=str(resource_review),
        primary_candidate_admission=False,scientific_energies_transferred=False,
        large_polymer_or_Hessian_resource_qualification=False))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('construction','resource-review','out'):p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--selection',action='append',required=True);a=p.parse_args()
    from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity
    with allocation_guard():host_capacity(1);os.sched_setaffinity(0,cpu_allocation(1))
    prepare(a.construction.resolve(),a.resource_review.resolve(),a.selection,a.out.resolve())


if __name__=='__main__':main()
