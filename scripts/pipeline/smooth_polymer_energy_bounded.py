#!/usr/bin/env python
"""Use a measured hybrid-energy RAM class after the preserved 32-GiB pilot.

Only the host-memory limits change. The original alternate-model energy driver,
coordinates, electronic equations, checkpointing and physical gates are retained.
"""
import argparse
import inspect
import json
import math
from pathlib import Path

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api

ORIGINAL = ROOT / 'scripts/pipeline/smooth_polymer_energy.py'


def prepare(geometry, review, resource, failed, out):
    original = api.load_module('smooth_energy_preserved_first_driver', ORIGINAL)
    source = api.read_plan(geometry)
    point = json.loads(resource.read_text())
    if point['status'] != 'electronic_resource_point_independently_reviewed':
        raise ValueError('The larger RAM bound requires an independently reviewed hybrid-energy point')
    api.verify_sources(point['source_sha256'])
    reference = Path(point['source']); reference_plan = api.read_plan(reference)
    task, other = source['tasks'][0], reference_plan['tasks'][0]
    classifier = api.load_module('smooth_energy_RAM_graph_class', ROOT / 'scripts/pipeline/prepare_conformer_geometry_batch.py')
    if (any(task[k] != other[k] for k in ('atoms', 'active_basis_functions', 'charge', 'mult'))
            or classifier.constitution(task['smiles']) != classifier.constitution(other['smiles'])):
        raise ValueError('The measured electronic RAM class has another constitution/state/basis size')
    raw_point = reference / 'energy/pilot/energy_and_gradient.json'
    electronic = json.loads(raw_point.read_text())['result']
    from redox.core.protocol import ACTIVE_SP
    if electronic['energy_protocol'] != ACTIVE_SP or not electronic['converged_smd'] or not electronic['converged_gas']:
        raise ValueError('The measured RAM point has a different hybrid single-point protocol')
    bound = point['host_RSS_bound_GiB']
    if (not math.isfinite(bound) or not 32 < bound <= 128
            or not 32*2**30 < point['peak_host_RSS_bytes'] <= bound*2**30
            or reference_plan['execution_workspace_GiB'] != 28):
        raise ValueError('The observed hybrid fitting memory does not support this new bound')
    failed_status = failed / 'experimental_energy/status.json'
    failure = json.loads(failed_status.read_text())
    guard_path = failed / 'host_guard/status.json'; guard = json.loads(guard_path.read_text())
    owner = api.load_module('smooth_energy_restart_stopped_owner', ROOT / 'scripts/pipeline/run_polymer_spacing.py')
    if (failure['status'] != 'experimental_hybrid_energy_failed_preserved'
            or owner.owner_alive(failure['owner']) or guard['maximum_RSS_GiB'] != 32
            or guard['peak_RSS_bytes'] <= 32*2**30 or guard['free_RAM_bytes'] < 32*2**30
            or api.digest(failed / 'input.xyz') != api.digest(Path(json.loads(review.read_text())['state_dir']) / 'opt.xyz')):
        raise ValueError('The preserved attempt is not this geometry with an exceeded 32-GiB RAM guard')
    original.prepare(geometry, review, out)
    plan = api.read_plan(out)
    for p in (resource, reference / 'plan.json', raw_point, failed / 'plan.json', failed_status, guard_path):
        plan['input_sha256'][str(p)] = api.digest(p)
    plan['input_sha256'].update(point['source_sha256'])
    plan['producer_sha256'][str(Path(__file__).resolve())] = api.digest(__file__)
    plan.update(execution_driver=str(Path(__file__).resolve()), maximum_host_RSS_GiB=bound,
        electronic_RAM_class_review=str(resource), preserved_failed_attempt=str(failed),
        resource_scope='Measured same-constitution legacy hybrid RAM bound only; alternate-model pilot still prospective')
    atomic_write_json(out / 'plan.json', plan)
    print(json.dumps(dict(status='measured_RAM_class_alternate_energy_restart_prepared',
        host_RSS_bound_GiB=bound, resource_point_geometry_or_energy_transferred=False)))


def run(campaign):
    module = api.load_module('smooth_energy_unchanged_equations', ORIGINAL)
    body = inspect.getsource(module.run)
    replacements = {
        'host_capacity(2, minimum_ram_gib=64)': "host_capacity(2, minimum_ram_gib=plan['maximum_host_RSS_GiB']+32)",
        "guard['maximum_RSS_GiB'] != 32": "guard['maximum_RSS_GiB'] != plan['maximum_host_RSS_GiB']"}
    for old, new in replacements.items():
        if body.count(old) != 1: raise ValueError('The preserved energy/RAM guard template changed')
        body = body.replace(old, new, 1)
    scope = dict(vars(module)); local = {}
    exec(compile(body, str(Path(__file__).resolve()), 'exec'), scope, local)
    local['run'](campaign)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('prepare', 'run'))
    parser.add_argument('--campaign', type=Path, required=True)
    for name in ('geometry', 'review', 'resource', 'failed-attempt'): parser.add_argument('--'+name, type=Path)
    args = parser.parse_args()
    if args.action == 'prepare':
        if any(getattr(args, k) is None for k in ('geometry', 'review', 'resource', 'failed_attempt')):
            parser.error('prepare requires geometry/review/resource/failed-attempt')
        prepare(args.geometry.resolve(), args.review.resolve(), args.resource.resolve(),
                args.failed_attempt.resolve(), args.campaign.resolve())
    else: run(args.campaign.resolve())


if __name__ == '__main__': main()
