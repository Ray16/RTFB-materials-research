#!/usr/bin/env python
"""Prepare separate conformer geometries from a measured graph resource class.

Only execution resources are shared across starts/configurations. Coordinates,
chemical identities, optimizer results and all physical checks remain individual.
"""
import argparse
import json
import os
from pathlib import Path
import shutil
import time

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api

GEOMETRY = ROOT / 'scripts/pipeline/optimize_polymer_geometry.py'
GC = ROOT / 'scripts/pipeline/smd_gc_workspace.py'


def constitution(smiles):
    from rdkit import Chem
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise ValueError('Invalid molecular graph')
    Chem.RemoveStereochemistry(molecule)
    return Chem.MolToSmiles(molecule)


def resource_compatible(task, reference, protocol, reference_protocol):
    if (protocol != reference_protocol or (task['charge'], task['mult']) != (0, 1)
            or any(task.get(k) != reference.get(k) for k in ('atoms', 'charge', 'mult', 'active_basis_functions'))
            or constitution(task['smiles']) != constitution(reference['smiles'])):
        raise ValueError('Conformer geometry resource class differs in graph, state, basis size or protocol')


def prepare(sources, reference, destination, exclude):
    parent = api.read_plan(reference)
    api.verify_sources(parent['producer_sha256']); api.verify_sources(parent['input_sha256'])
    report_path, guard_path = reference / 'geometry/status.json', reference / 'host_guard/status.json'
    report, guard = (json.loads(p.read_text()) for p in (report_path, guard_path))
    if (report['status'] != 'stationary_geometry_checked_active_parent_and_native_minimum_pending'
            or guard['status'] != 'program_returned_review_required'
            or not 0 < guard['peak_RSS_bytes'] <= 4 * 2**30
            or api.load_module('batch_reference_owner', ROOT / 'scripts/pipeline/run_polymer_spacing.py').owner_alive(report['owner'])):
        raise ValueError('A finished, physically checked geometry with measured RSS below 4 GiB is required')
    api.verify_sources(report['output_sha256'])
    if len(parent['tasks']) != 1:
        raise ValueError('One measured geometry resource class is required')
    declaration = []
    destination.mkdir(parents=True, exist_ok=False)
    for source in sources:
        old = api.read_plan(source)
        api.verify_sources(old['producer_sha256']); api.verify_sources(old['input_sha256'])
        for task in old['tasks']:
            if task['key'] in exclude:
                continue
            resource_compatible(task, parent['tasks'][0], old['quantum_settings'], parent['quantum_settings'])
            out = destination / task['key']
            out.mkdir()
            shutil.copytree(source / 'controls', out / 'controls')
            ready_dir = out / 'live_inputs' / task['key']
            shutil.copytree(source / 'live_inputs' / task['key'], ready_dir)
            ready_path, xyz = ready_dir / 'task.json', ready_dir / 'input.xyz'
            ready = json.loads(ready_path.read_text()); ready['input_xyz'] = str(xyz)
            ready['input_sha256'][str(xyz)] = api.digest(xyz)
            atomic_write_json(ready_path, ready)
            for stage in ('energy', 'thermal'):
                for name in ('tasks', 'done', 'claims', 'failures', 'workers', 'pilot', 'launches'):
                    (out / stage / name).mkdir(parents=True)
            inputs = dict(old['input_sha256'])
            for p in (source / 'plan.json', reference / 'plan.json', report_path, guard_path, ready_path, xyz):
                inputs[str(p)] = api.digest(p)
            inputs.update(report['output_sha256'])
            producers = dict(old['producer_sha256'])
            producers.update(parent['producer_sha256'])
            for p in (Path(__file__).resolve(), GEOMETRY, GC):
                producers[str(p)] = api.digest(p)
            plan = dict(old, tasks=[task], input_sha256=inputs, producer_sha256=producers,
                raw_root=str(ROOT / 'calcs/polymer_models' / destination.name / task['key']),
                geometry_resource_reference=str(reference), maximum_host_RSS_GiB=16,
                minimum_free_RAM_GiB=32,
                geometry_resource_policy=dict(minimum_GPU_MiB=32768, maximum_peak_GPU_MiB=31744,
                    minimum_host_RAM_GiB=48),
                geometry_only_execution=True, energy_only_execution=True, native_pipeline_enabled=False,
                externally_managed_stages=True, prepared=time.time(),
                purpose='Separate same-constitution conformer optimization; measured resource class only, each physical gate and peak enforced anew')
            atomic_write_json(out / 'plan.json', plan)
            declaration.append(dict(key=task['key'], campaign=str(out), source=str(source),
                plan_sha256=api.digest(out / 'plan.json'), status='prepared_geometry_only'))
    manifest = dict(status='conformer_geometry_batch_prepared', entries=declaration,
        measured_reference=str(reference), measured_reference_peak_RSS_bytes=guard['peak_RSS_bytes'],
        maximum_RSS_GiB=16, minimum_free_RAM_GiB=32, minimum_GPU_MiB=32768,
        worker=str(Path(__file__).resolve()), quantum_results_transferred=False,
        source_sha256={str(p): api.digest(p) for p in (Path(__file__).resolve(), GEOMETRY, GC)})
    atomic_write_json(destination / 'batch.json', manifest)
    print(json.dumps(dict(prepared=len(declaration), maximum_RSS_GiB=16)))


def run(campaign):
    from redox.workflows.resources import host_capacity
    plan = api.read_plan(campaign)
    api.verify_sources(plan['producer_sha256']); api.verify_sources(plan['input_sha256'])
    guard = json.loads((campaign / 'host_guard/status.json').read_text())
    if (guard['pid'] != os.getpid() or guard['maximum_RSS_GiB'] != 16 or guard['minimum_free_GiB'] != 32):
        raise ValueError('The declared same-PID host guard must be active')
    geometry = api.load_module('same_graph_conformer_geometry_physical_gates', GEOMETRY)
    # The retained driver has a 192-GiB floor from its older composite resource
    # profile. This geometry-only class uses measured RSS plus the guarded floor.
    geometry.host_capacity = lambda threads, **kwargs: host_capacity(threads, minimum_ram_gib=48)
    gc = api.load_module('same_graph_conformer_geometry_GC', GC)
    metrics = dict(collections=0, unreachable_objects_collected=0, started=time.time())
    try:
        with gc.collect_before_build(metrics):
            geometry.run(campaign)
    finally:
        metrics['finished'] = time.time()
        atomic_write_json(campaign / 'geometry_GC.json', metrics)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('prepare', 'run'))
    parser.add_argument('--source', type=Path, action='append')
    parser.add_argument('--reference', type=Path)
    parser.add_argument('--out', type=Path)
    parser.add_argument('--campaign', type=Path)
    parser.add_argument('--exclude-key', action='append', default=[])
    args = parser.parse_args()
    if args.action == 'run':
        if not args.campaign:
            parser.error('run requires --campaign')
        run(args.campaign.resolve())
    else:
        if not args.source or not args.reference or not args.out:
            parser.error('prepare requires --source, --reference and --out')
        from redox.workflows.resources import allocation_guard, cpu_allocation, host_capacity
        with allocation_guard():
            host_capacity(1); os.sched_setaffinity(0, cpu_allocation(1))
        prepare([p.resolve() for p in args.source], args.reference.resolve(), args.out.resolve(), set(args.exclude_key))


if __name__ == '__main__':
    main()
