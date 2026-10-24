#!/usr/bin/env python
"""Separate geometry from hybrid energies after an exact-input resource point.

The point establishes bounded execution at the input, not optimizer convergence.
The complete geometry stage still enforces topology, conservative force checks
and its observed GPU bound; native Hessians and active energies remain separate.
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
CHECKER = ROOT / 'scripts/validation/qualify_gc_workspace.py'


def point_evidence(pilot, review):
    plan = api.read_plan(pilot)
    api.verify_sources(plan['producer_sha256'])
    api.verify_sources(plan['input_sha256'])
    proof = json.loads((review / 'summary.json').read_text())
    if (proof['status'] != 'electronic_resource_point_independently_reviewed'
            or proof['source'] != str(pilot) or proof['G_computed']
            or proof['optimized_geometry_qualified']):
        raise ValueError('An independently reviewed exact-input resource point is required')
    api.verify_sources(proof['source_sha256'])
    if len(plan['tasks']) != 1:
        raise ValueError('One exact graph and conformer is required')
    task = plan['tasks'][0]
    if (task['charge'], task['mult']) != (0, 1) or proof['task_key'] != task['key']:
        raise ValueError('The point and neutral-singlet task must match')
    point_path = pilot / 'energy/pilot/energy_and_gradient.json'
    point = json.loads(point_path.read_text())
    xyz = Path(point['result']['xyz'])
    original = pilot / 'live_inputs' / task['key'] / 'input.xyz'
    if api.digest(xyz) != api.digest(original):
        raise ValueError('Resource point must use the exact staged input coordinates')
    from redox.core.protocol import ACTIVE_SP, geom_sha1
    if (point['result']['energy_geom_sha1'] != geom_sha1(xyz)
            or point['result']['energy_protocol'] != ACTIVE_SP):
        raise ValueError('The resource point geometry or active-energy protocol differs')
    return plan, proof, original


def prepare(pilot, review, controls, campaign):
    old, proof, original = point_evidence(pilot, review)
    evidence = api.load_module('point_geometry_GC_controls', CHECKER).checked_controls(controls)
    task = old['tasks'][0]
    campaign.mkdir(parents=True, exist_ok=False)
    shutil.copytree(pilot / 'controls', campaign / 'controls')
    shutil.copytree(pilot / 'live_inputs', campaign / 'live_inputs')
    for stage in ('energy', 'thermal'):
        for name in ('tasks', 'done', 'claims', 'failures', 'workers', 'pilot', 'launches'):
            (campaign / stage / name).mkdir(parents=True)
    ready_path = campaign / 'live_inputs' / task['key'] / 'task.json'
    ready = json.loads(ready_path.read_text())
    xyz = ready_path.parent / 'input.xyz'
    ready['input_xyz'] = str(xyz)
    ready['input_sha256'][str(xyz)] = api.digest(xyz)
    atomic_write_json(ready_path, ready)
    inputs = dict(old['input_sha256'])
    inputs.update(evidence)
    for path in (pilot / 'plan.json', review / 'summary.json', original, ready_path, xyz):
        inputs[str(path)] = api.digest(path)
    sources = dict(old['producer_sha256'])
    sources.update(evidence)
    for path in (Path(__file__).resolve(), GEOMETRY, GC, CHECKER):
        sources[str(path)] = api.digest(path)
    plan = dict(old, raw_root=str(ROOT / 'calcs/polymer_models' / campaign.parent.name / campaign.name),
        input_sha256=inputs, producer_sha256=sources, geometry_resource_pilot=str(pilot),
        geometry_resource_review=str(review / 'summary.json'),
        geometry_resource_policy=dict(minimum_GPU_MiB=32768, maximum_peak_GPU_MiB=31744,
            minimum_host_RAM_GiB=192),
        maximum_host_RSS_GiB=proof['host_RSS_bound_GiB'], minimum_free_RAM_GiB=32,
        geometry_only_execution=True, energy_only_execution=True, native_pipeline_enabled=False,
        externally_managed_stages=True, prepared=time.time(),
        purpose='Exact input-qualified geometry optimization; independent force/graph checks required; active energy/native G are separate')
    atomic_write_json(campaign / 'plan.json', plan)


def run(campaign):
    plan = api.read_plan(campaign)
    api.verify_sources(plan['producer_sha256'])
    api.verify_sources(plan['input_sha256'])
    guard = json.loads((campaign / 'host_guard/status.json').read_text())
    if (guard['pid'] != os.getpid() or guard['maximum_RSS_GiB'] != plan['maximum_host_RSS_GiB']
            or guard['minimum_free_GiB'] != plan['minimum_free_RAM_GiB']):
        raise ValueError('The declared owned host-resource guard must be active')
    module = api.load_module('point_qualified_geometry_physical_pipeline', GEOMETRY)
    gc = api.load_module('point_qualified_geometry_GC', GC)
    report = dict(collections=0, unreachable_objects_collected=0, started=time.time(),
        scientific_model_changed=False, source_sha256={str(p): api.digest(p) for p in (GC, Path(__file__).resolve())})
    try:
        with gc.collect_before_build(report):
            module.run(campaign)
    finally:
        report['finished'] = time.time()
        atomic_write_json(campaign / 'geometry_GC.json', report)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('prepare', 'run'))
    for name in ('pilot', 'review', 'controls', 'campaign'):
        parser.add_argument('--' + name, type=Path)
    args = parser.parse_args()
    if args.campaign is None:
        parser.error('--campaign is required')
    if args.action == 'prepare':
        if any(getattr(args, name) is None for name in ('pilot', 'review', 'controls')):
            parser.error('prepare requires --pilot, --review and --controls')
        from redox.workflows.resources import allocation_guard, cpu_allocation, host_capacity
        with allocation_guard():
            host_capacity(1)
            os.sched_setaffinity(0, cpu_allocation(1))
        prepare(args.pilot.resolve(), args.review.resolve(), args.controls.resolve(), args.campaign.resolve())
    else:
        run(args.campaign.resolve())


if __name__ == '__main__':
    main()
