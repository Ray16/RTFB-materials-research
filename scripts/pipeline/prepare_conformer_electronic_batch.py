#!/usr/bin/env python
"""Bind canonical hybrid energies to checked conformer geometries as they finish.

The measured point supplies an execution resource class, never an energy or
geometry. Every conformer receives its own normal gas/SMD single-point record.
"""
import argparse
import json
from pathlib import Path
import shutil
import time

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api

DRIVER = ROOT / 'scripts/pipeline/finish_polymer_electronic_gpu.py'
GC = ROOT / 'scripts/pipeline/smd_gc_workspace.py'
CHECKER = ROOT / 'scripts/validation/qualify_gc_workspace.py'
CLASSIFIER = ROOT / 'scripts/pipeline/prepare_conformer_geometry_batch.py'


def declare(geometry_batch, pilot, review, controls, out):
    upstream = json.loads((geometry_batch / 'batch.json').read_text())
    api.verify_sources(upstream['source_sha256'])
    point = api.read_plan(pilot)
    api.verify_sources(point['producer_sha256']); api.verify_sources(point['input_sha256'])
    proof_path = review / 'summary.json'; proof = json.loads(proof_path.read_text())
    if proof['status'] != 'electronic_resource_point_independently_reviewed' or proof['source'] != str(pilot):
        raise ValueError('The measured resource class needs independent review')
    api.verify_sources(proof['source_sha256'])
    witness = api.load_module('conformer_energy_GC_controls', CHECKER).checked_controls(controls)
    classifier = api.load_module('conformer_energy_graph_resource_class', CLASSIFIER)
    entries = []
    for entry in upstream['entries']:
        geometry = Path(entry['campaign']); plan = api.read_plan(geometry)
        if api.digest(geometry / 'plan.json') != entry['plan_sha256'] or len(plan['tasks']) != 1:
            raise ValueError('The declared exact conformer geometry changed')
        classifier.resource_compatible(plan['tasks'][0], point['tasks'][0],
            plan['quantum_settings'], point['quantum_settings'])
        entries.append(dict(key=entry['key'], campaign=str(out / entry['key']),
            geometry_campaign=str(geometry), geometry_plan_sha256=entry['plan_sha256'],
            plan_sha256=None))
    out.mkdir(parents=True, exist_ok=False)
    sources = {**witness, str(pilot / 'plan.json'): api.digest(pilot / 'plan.json'),
        str(proof_path): api.digest(proof_path), str(geometry_batch / 'batch.json'): api.digest(geometry_batch / 'batch.json')}
    sources.update({str(p): api.digest(p) for p in (Path(__file__).resolve(), DRIVER, GC, CHECKER, CLASSIFIER)})
    manifest = dict(status='deferred_conformer_electronics_declared', stage='electronic',
        entries=entries, worker=str(DRIVER), preparer=str(Path(__file__).resolve()),
        maximum_RSS_GiB=proof['host_RSS_bound_GiB'], minimum_free_RAM_GiB=32, minimum_GPU_MiB=32768,
        pilot=str(pilot), resource_review=str(proof_path), controls=str(controls),
        source_sha256=sources, created=time.time(), quantum_results_transferred=False)
    atomic_write_json(out / 'batch.json', manifest)
    print(json.dumps(dict(declared=len(entries), stage='electronic', geometry_dependency=True)))


def prepare_ready(entry, manifest):
    geometry, out = Path(entry['geometry_campaign']), Path(entry['campaign'])
    if api.digest(geometry / 'plan.json') != entry['geometry_plan_sha256']:
        raise ValueError('Upstream geometry declaration changed')
    path = geometry / 'geometry/status.json'
    if not path.exists():
        return 'upstream_pending'
    parent = json.loads(path.read_text())
    if 'finished' not in parent:
        return 'upstream_pending'
    if parent['status'] != 'stationary_geometry_checked_active_parent_and_native_minimum_pending':
        return 'upstream_failed_preserved'
    alive = api.load_module('conformer_energy_upstream_owner', ROOT / 'scripts/pipeline/run_polymer_spacing.py').owner_alive
    if alive(parent['owner']):
        return 'upstream_pending'
    if (out / 'plan.json').exists():
        validate_prepared(entry, manifest)
        return 'prepared'
    original = api.read_plan(geometry); point = api.read_plan(Path(manifest['pilot']))
    proof = json.loads(Path(manifest['resource_review']).read_text())
    for p in (original, point):
        api.verify_sources(p['producer_sha256']); api.verify_sources(p['input_sha256'])
    api.verify_sources(parent['output_sha256']); api.verify_sources(proof['source_sha256'])
    classifier = api.load_module('checked_energy_conformer_resource_class', CLASSIFIER)
    classifier.resource_compatible(original['tasks'][0], point['tasks'][0], original['quantum_settings'], point['quantum_settings'])
    witness = api.load_module('checked_energy_conformer_GC_controls', CHECKER).checked_controls(Path(manifest['controls']))
    merge = api.load_module('checked_energy_conformer_hash_join', DRIVER).merge_hashes
    inputs = merge(original['input_sha256'], point['input_sha256'], parent['output_sha256'], proof['source_sha256'], witness)
    for p in (geometry / 'plan.json', Path(manifest['pilot']) / 'plan.json', path, Path(manifest['resource_review'])):
        inputs[str(p)] = api.digest(p)
    sources = merge(original['producer_sha256'], point['producer_sha256'], witness)
    sources.update({str(p): api.digest(p) for p in (Path(__file__).resolve(), DRIVER, GC, CHECKER,
        CLASSIFIER, ROOT / 'src/redox/workflows/split_electronic.py')})
    plan = dict(original, raw_root=str(ROOT / 'calcs/polymer_models' / out.parent.name / out.name),
        producer_sha256=sources, input_sha256=inputs, source_geometry_report=parent,
        GPU_resource_qualification=manifest['pilot'], GPU_resource_review=manifest['resource_review'],
        resource_qualification_scope='Same constitution/state/basis execution class; no coordinates or energies transferred',
        maximum_RSS_GiB=manifest['maximum_RSS_GiB'], minimum_free_RAM_GiB=32,
        maximum_peak_GPU_MiB=point['maximum_memory_pilot_MiB'], execution_workspace_GiB=point['execution_workspace_GiB'],
        energy_only_execution=True, geometry_only_execution=False, native_pipeline_enabled=False,
        externally_managed_stages=True, created=time.time(),
        purpose='Fresh canonical hybrid gas/SMD SP on exact independently stationary conformer coordinates')
    out.mkdir(parents=True, exist_ok=False)
    shutil.copytree(geometry / 'controls', out / 'controls')
    shutil.copytree(geometry / 'live_inputs', out / 'live_inputs')
    key = original['tasks'][0]['key']; ready_path = out / 'live_inputs' / key / 'task.json'
    ready = json.loads(ready_path.read_text()); xyz = ready_path.parent / 'input.xyz'
    ready['input_xyz'] = str(xyz); ready['input_sha256'][str(xyz)] = api.digest(xyz)
    atomic_write_json(ready_path, ready)
    inputs[str(ready_path)] = api.digest(ready_path); inputs[str(xyz)] = api.digest(xyz)
    for stage in ('energy', 'thermal'):
        for name in ('tasks', 'done', 'claims', 'failures', 'workers', 'pilot', 'launches'):
            (out / stage / name).mkdir(parents=True)
    atomic_write_json(out / 'plan.json', plan)
    validate_prepared(entry, manifest)
    return 'prepared'


def validate_prepared(entry, manifest):
    geometry, out = Path(entry['geometry_campaign']), Path(entry['campaign'])
    if api.digest(geometry / 'plan.json') != entry['geometry_plan_sha256']:
        raise ValueError('Upstream geometry declaration changed')
    original, plan = api.read_plan(geometry), api.read_plan(out)
    parent_path = geometry / 'geometry/status.json'; parent = json.loads(parent_path.read_text())
    if (parent['status'] != 'stationary_geometry_checked_active_parent_and_native_minimum_pending'
            or plan['tasks'] != original['tasks'] or plan['source_geometry_report'] != parent
            or plan['GPU_resource_qualification'] != manifest['pilot']
            or plan['input_sha256'].get(str(parent_path)) != api.digest(parent_path)
            or plan['maximum_RSS_GiB'] != manifest['maximum_RSS_GiB']):
        raise ValueError('The electronic plan differs from its declared exact stationary parent')
    api.verify_sources(plan['producer_sha256']); api.verify_sources(plan['input_sha256'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('geometry-batch', 'pilot', 'review', 'controls', 'out'):
        parser.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args()
    declare(args.geometry_batch.resolve(), args.pilot.resolve(), args.review.resolve(), args.controls.resolve(), args.out.resolve())


if __name__ == '__main__':
    main()
