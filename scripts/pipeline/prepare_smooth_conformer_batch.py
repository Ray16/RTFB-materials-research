#!/usr/bin/env python
"""Prepare same-constitution smooth-model geometries after a checked pilot.

Only measured resource bounds transfer. Coordinates, electronic calculations,
force checks and subsequent thermochemistry remain independent for every start.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import time

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api


def prepare(sources,reference,review,out):
    from ase.io import read
    model_path=ROOT/'scripts/pipeline/smooth_cds_control.py'
    worker=ROOT/'scripts/pipeline/run_smooth_geometry.py'
    coordinator=ROOT/'scripts/pipeline/coordinate_smooth_geometry.py'
    model=api.load_module('smooth_batch_model_definition',model_path)
    classifier=api.load_module('smooth_batch_resource_identity',ROOT/'scripts/pipeline/prepare_conformer_geometry_batch.py')
    checked=json.loads(review.read_text())
    if checked['status']!='alternate_smooth_geometry_and_resource_pilot_independently_checked':
        raise ValueError('Independent alternate-model resource/geometry review required')
    api.verify_sources(checked['source_sha256'])
    baseline=api.read_plan(reference);task_reference=baseline['tasks'][0]
    if checked['campaign']!=str(reference) or checked['method']!=model.MODEL:
        raise ValueError('Measured reference has another campaign or model')
    identity=hashlib.sha256(json.dumps(model.MODEL,sort_keys=True).encode()).hexdigest()
    if checked['portable_method_id']!=identity or checked['peak_host_RSS_bytes']>4*2**30:
        raise ValueError('Measured model/class bound is not supported')
    out.mkdir(parents=True,exist_ok=False);entries=[]
    for source in sources:
        original=api.read_plan(source);api.verify_sources(original['producer_sha256']);api.verify_sources(original['input_sha256'])
        for task in original['tasks']:
            if ((task['charge'],task['mult'])!=(0,1)
                    or any(task[k]!=task_reference[k] for k in ('atoms','active_basis_functions'))
                    or classifier.constitution(task['smiles'])!=classifier.constitution(task_reference['smiles'])):
                raise ValueError('A different constitution/state/basis needs its own new-model resource pilot')
            case=out/task['key'];case.mkdir()
            ready=case/'live_inputs'/task['key'];shutil.copytree(source/'live_inputs'/task['key'],ready)
            xyz=ready/'input.xyz';record=json.loads((ready/'task.json').read_text())
            record['input_xyz']=str(xyz);record['input_sha256'][str(xyz)]=api.digest(xyz)
            atomic_write_json(ready/'task.json',record)
            for stage in ('energy','thermal'):
                for name in ('tasks','done','claims','failures','workers','pilot','launches'):
                    (case/stage/name).mkdir(parents=True)
            inputs=dict(original['input_sha256']);inputs.update({str(p):api.digest(p) for p in
                (source/'plan.json',reference/'plan.json',review,xyz,ready/'task.json',ROOT/'environments/cutensor22/runtime.json')})
            inputs.update(checked['source_sha256'])
            producers=dict(original['producer_sha256']);producers.update(checked['source_sha256'])
            for p in (Path(__file__).resolve(),worker,coordinator,model_path):producers[str(p)]=api.digest(p)
            symbols=read(xyz).get_chemical_symbols();axes=[[0,1]]
            hetero=next((i for i,s in enumerate(symbols) if s not in ('C','H')),None)
            if hetero is not None:axes.append([hetero,2])
            plan=dict(original,tasks=[task],input_sha256=inputs,producer_sha256=producers,
                execution_driver=str(worker),runtime_profile=str(ROOT/'environments/cutensor22/runtime.json'),
                raw_root=str(ROOT/'calcs/polymer_models'/out.name/task['key']),maximum_host_RSS_GiB=32,
                quantum_settings=dict(original['quantum_settings'],geometry_SCF_conv_tol=1e-12),
                geometry_resource_policy=dict(minimum_GPU_MiB=32768,maximum_peak_GPU_MiB=31744,minimum_host_RAM_GiB=64),
                experimental_model=model.MODEL,experimental_method_id=identity,
                full_quantum_controls=baseline['full_quantum_controls'],extra_force_axes=axes,
                primary_candidate_admission=False,native_pipeline_enabled=False,geometry_only_execution=True,
                energy_only_execution=True,externally_managed_stages=True,
                measured_alternate_model_reference=str(reference),resource_review=str(review),created=time.time(),
                purpose='Independent alternate-model conformer geometry; resources only inherited, no energies/G copied')
            atomic_write_json(case/'plan.json',plan)
            entries.append(dict(key=task['key'],campaign=str(case),plan_sha256=api.digest(case/'plan.json')))
    atomic_write_json(out/'batch.json',dict(status='alternate_model_conformer_geometries_prepared',stage='geometry',
        entries=entries,worker=str(worker),maximum_RSS_GiB=32,minimum_free_RAM_GiB=32,minimum_GPU_MiB=32768,
        experimental_model=model.MODEL,primary_candidate_admission=False,
        source_sha256={str(p):api.digest(p) for p in (Path(__file__).resolve(),worker,coordinator,model_path,
            ROOT/'scripts/pipeline/coordinate_geometry_batch.py',ROOT/'environments/cutensor22/runtime.json')},
        quantum_results_transferred=False))
    print(json.dumps(dict(prepared=len(entries),primary_candidate_admission=False)))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',type=Path,action='append',required=True)
    for name in ('reference','review','out'):parser.add_argument('--'+name,type=Path,required=True)
    args=parser.parse_args();prepare([p.resolve() for p in args.source],args.reference.resolve(),args.review.resolve(),args.out.resolve())


if __name__=='__main__':main()
