#!/usr/bin/env python
"""Bind fresh alternate-model energies to independently reviewed geometries.

Resource evidence transfers only within the measured constitution/state/basis
and model class. A completed energy is reused only for its identical declared
geometry; it is never copied onto another conformer.
"""
import argparse
import json
import math
from pathlib import Path
import time

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api

WORKER = ROOT / 'scripts/pipeline/smooth_polymer_energy_bounded.py'
COORDINATOR = ROOT / 'scripts/pipeline/coordinate_smooth_electronic_batch.py'


def resource_compatible(plan, point, review):
    classifier = api.load_module('alternate_electronic_constitution', ROOT / 'scripts/pipeline/prepare_conformer_geometry_batch.py')
    task, other = plan['tasks'][0], point['tasks'][0]
    if (len(plan['tasks']) != 1 or len(point['tasks']) != 1 or (task['charge'], task['mult']) != (0, 1)
            or any(task[k] != other[k] for k in ('atoms', 'active_basis_functions', 'charge', 'mult'))
            or classifier.constitution(task['smiles']) != classifier.constitution(other['smiles'])
            or plan['experimental_method_id'] != point['experimental_method_id']
            or plan['experimental_model'] != point['experimental_model']
            or plan['quantum_settings'] != point['quantum_settings']
            or plan['runtime_profile'] != point['runtime_profile']):
        raise ValueError('The alternate electronic execution class differs in graph/state/basis/model/settings')
    bound = point['maximum_host_RSS_GiB']
    if (review['status'] != 'alternate_electronic_geometry_and_resource_checked_native_G_pending'
            or review['portable_method_id'] != point['experimental_method_id']
            or review['native_Hessian_or_G_qualified'] or review['primary_candidate_admission']
            or not math.isfinite(bound) or not 0 < review['peak_host_RSS_bytes'] <= bound*2**30
            or not 0 < review['observed_peak_GPU_MiB'] <= point['maximum_peak_GPU_MiB']):
        raise ValueError('Completed alternate electronic resource evidence is invalid')


def declare(geometry_batch, geometry_reviews, reference, resource_review, out):
    source_path = geometry_batch / 'batch.json'; upstream = json.loads(source_path.read_text())
    if upstream['stage'] != 'geometry' or upstream['primary_candidate_admission']:
        raise ValueError('Only a declared alternate-model geometry batch is supported')
    api.verify_sources(upstream['source_sha256'])
    point = api.read_plan(reference); review = json.loads(resource_review.read_text())
    api.verify_sources(point['producer_sha256']); api.verify_sources(point['input_sha256'])
    api.verify_sources(review['source_sha256'])
    if review['campaign'] != str(reference):
        raise ValueError('Resource review belongs to another energy point')
    entries = []
    for entry in upstream['entries']:
        geometry = Path(entry['campaign']); plan = api.read_plan(geometry)
        if api.digest(geometry / 'plan.json') != entry['plan_sha256']:
            raise ValueError('The declared source geometry plan changed')
        resource_compatible(plan, point, review)
        identical = str(geometry) == point['source_geometry_campaign']
        entries.append(dict(key=entry['key'], campaign=str(reference if identical else out / entry['key']),
            geometry_campaign=str(geometry), geometry_plan_sha256=entry['plan_sha256'],
            geometry_review=str(geometry_reviews / (geometry_batch.name+'__'+entry['key']+'.json')),
            existing_checked_reference=identical, plan_sha256=None))
    out.mkdir(parents=True, exist_ok=False)
    paths = [Path(__file__).resolve(), COORDINATOR, WORKER,
        ROOT/'scripts/pipeline/smooth_polymer_energy.py', ROOT/'scripts/pipeline/dispatch_smooth_energy_pilot.py',
        ROOT/'scripts/pipeline/coordinate_geometry_batch.py', ROOT/'scripts/validation/review_smooth_polymer_energy.py',
        ROOT/'environments/cutensor22/runtime.json', source_path, reference/'plan.json', resource_review]
    sources = {str(p):api.digest(p) for p in paths}
    atomic_write_json(out / 'batch.json', dict(status='deferred_alternate_electronic_batch_declared',
        stage='electronic', worker=str(WORKER), preparer=str(Path(__file__).resolve()), entries=entries,
        reference=str(reference), resource_review=str(resource_review),
        maximum_RSS_GiB=point['maximum_host_RSS_GiB'], minimum_free_RAM_GiB=32, minimum_GPU_MiB=32768,
        experimental_method_id=point['experimental_method_id'], primary_candidate_admission=False,
        native_G_computed=False, source_sha256=sources, created=time.time(),
        existing_reference_entries=sum(e['existing_checked_reference'] for e in entries),
        quantum_values_transferred_to_other_geometries=False))
    print(json.dumps(dict(declared=len(entries), identical_reference_reused=sum(e['existing_checked_reference'] for e in entries))))


def prepare_ready(entry, manifest):
    geometry, out = Path(entry['geometry_campaign']), Path(entry['campaign'])
    if api.digest(geometry / 'plan.json') != entry['geometry_plan_sha256']:
        raise ValueError('The declared geometry input changed')
    if entry['existing_checked_reference']:
        return 'existing_checked_reference'
    if (out / 'plan.json').exists():
        validate_prepared(entry, manifest); return 'prepared'
    if out.exists():
        return 'partial_preparation_preserved'
    status = geometry / 'geometry/status.json'
    if not status.exists(): return 'upstream_pending'
    state = json.loads(status.read_text())
    if 'finished' not in state: return 'upstream_pending'
    if state['status'] != 'experimental_smooth_geometry_checked_energy_and_native_QC_pending':
        return 'upstream_failed_preserved'
    owners = api.load_module('alternate_energy_geometry_ownership', ROOT / 'scripts/pipeline/run_polymer_spacing.py')
    if owners.owner_alive(state['owner']): return 'upstream_pending'
    review_path = Path(entry['geometry_review'])
    if not review_path.exists(): return 'upstream_pending'
    plan, point = api.read_plan(geometry), api.read_plan(Path(manifest['reference']))
    resource = json.loads(Path(manifest['resource_review']).read_text())
    api.verify_sources(resource['source_sha256']); resource_compatible(plan, point, resource)
    driver = api.load_module('alternate_energy_exact_geometry_preparation', ROOT / 'scripts/pipeline/smooth_polymer_energy.py')
    driver.prepare(geometry, review_path, out)
    fresh = api.read_plan(out)
    for p in (Path(manifest['reference']) / 'plan.json', Path(manifest['resource_review'])):
        fresh['input_sha256'][str(p)] = api.digest(p)
    fresh['input_sha256'].update(resource['source_sha256'])
    for p in (Path(__file__).resolve(), WORKER, COORDINATOR): fresh['producer_sha256'][str(p)] = api.digest(p)
    fresh.update(execution_driver=str(WORKER), maximum_host_RSS_GiB=manifest['maximum_RSS_GiB'],
        maximum_peak_GPU_MiB=point['maximum_peak_GPU_MiB'], prospective_resource_pilot=False,
        measured_alternate_model_electronic_reference=manifest['reference'],
        electronic_RAM_class_review=manifest['resource_review'],
        resource_scope='Measured same-constitution/state/basis/model class; no coordinates/energies transferred')
    atomic_write_json(out / 'plan.json', fresh)
    validate_prepared(entry, manifest)
    return 'prepared'


def validate_prepared(entry, manifest):
    geometry, out = Path(entry['geometry_campaign']), Path(entry['campaign'])
    if api.digest(geometry / 'plan.json') != entry['geometry_plan_sha256']:
        raise ValueError('Source geometry plan changed')
    plan, original = api.read_plan(out), api.read_plan(geometry)
    if (entry['existing_checked_reference'] or plan['tasks'] != original['tasks']
            or plan['source_geometry_campaign'] != str(geometry)
            or plan['execution_driver'] != str(WORKER)
            or plan['maximum_host_RSS_GiB'] != manifest['maximum_RSS_GiB']
            or plan['measured_alternate_model_electronic_reference'] != manifest['reference']
            or plan['experimental_method_id'] != manifest['experimental_method_id']
            or plan['primary_candidate_admission'] or plan['native_G_computed']):
        raise ValueError('Prepared alternate energy differs from its declared geometry/model/resource class')
    api.verify_sources(plan['producer_sha256']); api.verify_sources(plan['input_sha256'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('geometry-batch','geometry-reviews','reference','resource-review','out'):
        parser.add_argument('--'+name, type=Path, required=True)
    args = parser.parse_args()
    from redox.workflows.resources import allocation_guard, cpu_allocation, host_capacity
    import os
    with allocation_guard(): host_capacity(1); os.sched_setaffinity(0,cpu_allocation(1))
    declare(args.geometry_batch.resolve(), args.geometry_reviews.resolve(), args.reference.resolve(),
            args.resource_review.resolve(), args.out.resolve())


if __name__ == '__main__': main()
