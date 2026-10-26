#!/usr/bin/env python
"""Declare fresh native states after a same-class native resource pilot passes.

Electronic recovery options are tied to the same exact geometry. Only execution
bounds transfer between conformers; coordinates, energies and Hessians do not.
"""
import argparse
import json
import math
import os
from pathlib import Path
import time

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api

WORKER = ROOT / 'scripts/pipeline/smooth_polymer_native.py'


def resource_compatible(plan, point, review):
    classifier = api.load_module('smooth_native_constitution', ROOT / 'scripts/pipeline/prepare_conformer_geometry_batch.py')
    if len(plan['tasks']) != 1 or len(point['tasks']) != 1:
        raise ValueError('One exact state per native task is required')
    task, other = plan['tasks'][0], point['tasks'][0]
    if ((task['charge'], task['mult']) != (0, 1)
            or any(task[k] != other[k] for k in ('atoms', 'active_basis_functions', 'charge', 'mult'))
            or classifier.constitution(task['smiles']) != classifier.constitution(other['smiles'])
            or any(plan[k] != point[k] for k in ('experimental_method_id', 'experimental_model',
                                                'quantum_settings', 'runtime_profile'))):
        raise ValueError('Native resource class differs in graph/state/basis/model/settings/runtime')
    validate_review(review, point['experimental_method_id'])
    if (point['maximum_host_RSS_GiB'] != 80 or point['maximum_peak_GPU_MiB'] != 31744
            or not 0 < review['peak_host_RSS_bytes'] <= 80 * 2**30
            or not 0 < review['observed_peak_GPU_MiB'] <= 31744):
        raise ValueError('Same-class native resource measurements are missing or exceed their bounds')


def validate_review(review, method_id):
    if (review['status'] != 'alternate_native_state_independently_reconstructed'
            or review['experimental_method_id'] != method_id or review['primary_candidate_admission']
            or review['conformer_convergence_established']):
        raise ValueError('Independent native state has the wrong method or scientific scope')
    values = [review['G_eV'], review['minimum_frequency_cm'], review['maximum_sampled_column_error'],
              review['frequency_reconstruction_error_cm'], *review['G_by_cutoff_eV'].values(),
              *review['nuclear_reconstruction_errors_eV'].values()]
    if (set(review['G_by_cutoff_eV']) != {'25', '50', '100'}
            or set(review['nuclear_reconstruction_errors_eV']) != {'25', '50', '100'}
            or any(type(v) not in (int, float) or not math.isfinite(v) for v in values)
            or min(review['maximum_sampled_column_error'], review['frequency_reconstruction_error_cm'],
                   *review['nuclear_reconstruction_errors_eV'].values()) < 0
            or abs(review['G_eV'] - review['G_by_cutoff_eV']['50']) > 1e-8
            or review['minimum_frequency_cm'] < -5 or review['maximum_sampled_column_error'] > .01
            or review['frequency_reconstruction_error_cm'] > 1e-5
            or max(review['nuclear_reconstruction_errors_eV'].values()) > 1e-8):
        raise ValueError('Native free energies, modes or derivative reconstruction are invalid')


def declare(electronic_batches, reference, resource_review, out):
    point = api.read_plan(reference); witness = json.loads(resource_review.read_text())
    api.verify_sources(point['producer_sha256']); api.verify_sources(point['input_sha256'])
    api.verify_sources(witness['source_sha256'])
    if witness['campaign'] != str(reference):
        raise ValueError('Native resource review belongs to another state')
    resource_compatible(point, point, witness)
    groups = {}; sources = {}
    for batch in electronic_batches:
        path = batch / 'batch.json'; upstream = json.loads(path.read_text())
        if upstream['stage'] != 'electronic' or upstream['primary_candidate_admission']:
            raise ValueError('Only declared, separate-model electronic batches are supported')
        api.verify_sources(upstream['source_sha256']); sources[str(path)] = api.digest(path)
        for entry in upstream['entries']:
            geometry = Path(entry['geometry_campaign']); plan = api.read_plan(geometry)
            if api.digest(geometry / 'plan.json') != entry['geometry_plan_sha256']:
                raise ValueError('Electronic source geometry plan changed')
            resource_compatible(plan, point, witness)
            if entry['key'] in groups and groups[entry['key']]['geometry_campaign'] != str(geometry):
                raise ValueError('One conformer key refers to two different source geometries')
            identical = str(geometry) == point['source_geometry_campaign']
            row = groups.setdefault(entry['key'], dict(key=entry['key'],
                campaign=str(reference if identical else out / entry['key']), geometry_campaign=str(geometry),
                geometry_plan_sha256=entry['geometry_plan_sha256'], existing_checked_reference=identical,
                plan_sha256=None, electronic_options=[]))
            case = Path(entry['campaign'])
            review = Path(upstream['resource_review']) if entry['existing_checked_reference'] else case/'independent_energy_review.json'
            option = dict(campaign=str(case), review=str(review))
            if option not in row['electronic_options']: row['electronic_options'].append(option)
    if not groups or sum(r['existing_checked_reference'] for r in groups.values()) != 1:
        raise ValueError('One exact native reference must occur in the declared electronic cohorts')
    out.mkdir(parents=True, exist_ok=False)
    for path in (Path(__file__).resolve(), WORKER, ROOT/'scripts/pipeline/coordinate_smooth_native_batch.py',
                 ROOT/'scripts/validation/review_smooth_polymer_native.py',
                 ROOT/'scripts/pipeline/dispatch_smooth_energy_pilot.py',
                 ROOT/'scripts/pipeline/coordinate_geometry_batch.py',
                 ROOT/'scripts/fleet/preflight_cpu_dependencies.py', reference/'plan.json', resource_review,
                 ROOT/'environments/cutensor22/runtime.json'):
        sources[str(path)] = api.digest(path)
    manifest = dict(status='native_batch_declared_after_independent_same_class_pilot', stage='native',
        worker=str(WORKER), preparer=str(Path(__file__).resolve()), entries=list(groups.values()),
        reference=str(reference), resource_review=str(resource_review), maximum_RSS_GiB=80,
        minimum_free_RAM_GiB=32, minimum_GPU_MiB=32768,
        experimental_method_id=point['experimental_method_id'], source_sha256=sources,
        primary_candidate_admission=False, created=time.time(),
        thermochemistry_transferred_between_geometries=False, search_hits_are_weights=False)
    atomic_write_json(out/'batch.json', manifest)
    return manifest


def prepare_ready(entry, manifest):
    out = Path(entry['campaign']); geometry = Path(entry['geometry_campaign'])
    if api.digest(geometry/'plan.json') != entry['geometry_plan_sha256']:
        raise ValueError('Native source geometry plan changed')
    if entry['existing_checked_reference']: return 'existing_checked_reference'
    if (out/'plan.json').exists():
        validate_prepared(entry, manifest); return 'prepared'
    if out.exists(): return 'partial_preparation_preserved'
    owners = api.load_module('smooth_native_parent_ownership', ROOT/'scripts/pipeline/run_polymer_spacing.py')
    waiting = False
    for option in entry['electronic_options']:
        energy, review_path = Path(option['campaign']), Path(option['review'])
        status_path = energy/'experimental_energy/status.json'
        if not status_path.exists(): waiting = True; continue
        row = json.loads(status_path.read_text())
        if 'finished' not in row or owners.owner_alive(row['owner']): waiting = True; continue
        if row['status'] != 'experimental_hybrid_energy_finished_independent_review_pending': continue
        if not review_path.exists(): waiting = True; continue
        plan = api.read_plan(energy)
        if plan['source_geometry_campaign'] != str(geometry):
            raise ValueError('Electronic recovery belongs to a different geometry')
        point = api.read_plan(Path(manifest['reference']))
        resource = json.loads(Path(manifest['resource_review']).read_text())
        api.verify_sources(resource['source_sha256']); resource_compatible(plan, point, resource)
        driver = api.load_module('smooth_native_exact_parent_preparation', WORKER)
        driver.prepare(energy, review_path, out)
        fresh = api.read_plan(out)
        for path in (Path(manifest['reference'])/'plan.json', Path(manifest['resource_review'])):
            fresh['input_sha256'][str(path)] = api.digest(path)
        fresh['input_sha256'].update(resource['source_sha256'])
        fresh.update(prospective_resource_pilot=False, measured_native_reference=manifest['reference'],
            native_resource_review=manifest['resource_review'])
        atomic_write_json(out/'plan.json', fresh); validate_prepared(entry, manifest)
        return 'prepared'
    return 'upstream_pending' if waiting else 'upstream_failed_preserved'


def validate_prepared(entry, manifest):
    plan = api.read_plan(Path(entry['campaign']))
    api.verify_sources(plan['producer_sha256']); api.verify_sources(plan['input_sha256'])
    point = api.read_plan(Path(manifest['reference']))
    witness = json.loads(Path(manifest['resource_review']).read_text())
    resource_compatible(plan, point, witness)
    parents = {p['campaign'] for p in entry['electronic_options']}
    if (entry['existing_checked_reference'] or plan['source_electronic_case'] not in parents
            or plan['source_geometry_campaign'] != entry['geometry_campaign']
            or plan['execution_driver'] != str(WORKER) or plan['primary_candidate_admission']
            or plan['measured_native_reference'] != manifest['reference']
            or plan['experimental_native_settings'] != point['experimental_native_settings']):
        raise ValueError('Prepared native state differs from its declared exact parent/resource class')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--electronic-batch', type=Path, action='append', required=True)
    for name in ('reference', 'resource-review', 'out'): parser.add_argument('--'+name, type=Path, required=True)
    args = parser.parse_args()
    from redox.workflows.resources import allocation_guard, cpu_allocation, host_capacity
    with allocation_guard(): host_capacity(1); os.sched_setaffinity(0,cpu_allocation(1))
    manifest = declare([p.resolve() for p in args.electronic_batch], args.reference.resolve(),
                       args.resource_review.resolve(), args.out.resolve())
    print(json.dumps(dict(declared=len(manifest['entries']), primary_candidate_admission=False)))


if __name__ == '__main__': main()
