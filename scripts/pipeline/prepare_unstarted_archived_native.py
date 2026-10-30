#!/usr/bin/env python3
"""Stage an exact declared native state whose earlier wrapper never ran science."""
import argparse
import json
from pathlib import Path
import shutil

from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json
from redox.validation.archived_sources import ArchivedSources
from redox.workflows.archived_native_namespace import mount_command


def prepare(archive, prior, campaign, destination, inventory_path):
    inventory = json.loads(inventory_path.read_text())
    if len({r['host'] for r in inventory['nodes']}) != 9 or any(r['processes'] for r in inventory['nodes']):
        raise ValueError('Stopped-predictor evidence across the system-policy nodes required')
    reader = ArchivedSources(prior, archive, relative_to_prior=True)
    source = reader.path(campaign)
    plan = json.loads((source/'plan.json').read_text())
    task = plan['tasks'][0]
    if (len(plan['tasks']) != 1 or task['key'] not in inventory['keys']
            or (task['charge'], task['mult']) != (0, 1)
            or plan['experimental_model']['solvent'] != 'acetonitrile'
            or plan['execution_driver'] != str(prior/'scripts/pipeline/smooth_polymer_native.py')
            or plan['maximum_host_RSS_GiB'] != 80 or plan['maximum_peak_GPU_MiB'] != 31744):
        raise ValueError('Source is not the declared native neutral-singlet class')
    if (source/'experimental_native').exists() or reader.path(plan['raw_root']).exists():
        raise ValueError('This entry point cannot discard an existing native scientific attempt')
    review_path = plan['electronic_review']
    if review_path not in plan['input_sha256']:
        raise ValueError('Backing electronic review is not bound by the frozen plan')
    reader.verify_all(plan['producer_sha256'])
    reader.verify_all(plan['input_sha256'])
    review = json.loads(reader.path(review_path).read_text())
    if (review['status'] != 'alternate_electronic_geometry_and_resource_checked_native_G_pending'
            or review['native_Hessian_or_G_qualified']):
        raise ValueError('Independently checked electronic parent required')
    reader.verify_all(review['source_sha256'])
    if sha256_file(source/'input.xyz') != plan['scored_geometry_sha256']:
        raise ValueError('Native geometry differs from electronic parent')
    if destination.exists():
        raise FileExistsError('Preserve prior continuation destinations')
    attempt, raw = destination/'campaign', destination/'raw'
    attempt.mkdir(parents=True)
    raw.mkdir()
    for name in ('plan.json', 'input.xyz'):
        shutil.copyfile(source/name, attempt/name)
    spec = dict(status='prepared_native_state_after_stopped_pre_science_wrapper', archive=str(archive),
                prior_root=str(prior), original_campaign=str(campaign), original_raw_parent=str(Path(plan['raw_root']).parent),
                original_raw_root=plan['raw_root'], attempt=str(attempt), raw_parent=str(raw),
                plan_sha256=sha256_file(source/'plan.json'), geometry_sha256=plan['scored_geometry_sha256'],
                key=task['key'], method_id=plan['experimental_method_id'], maximum_host_RSS_GiB=80,
                minimum_free_GiB=32, historical_files_modified=False, scientific_settings_changed=False,
                inventory=str(inventory_path.resolve()), inventory_sha256=sha256_file(inventory_path),
                complete_native_G_supplied=False)
    mount_command(spec)
    reader.recheck()
    atomic_write_json(destination/'namespace.json', spec)
    return spec


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('archive', 'prior-root', 'campaign', 'destination', 'inventory'):
        p.add_argument('--'+name, type=Path, required=True)
    a = p.parse_args()
    result = prepare(a.archive.resolve(), a.prior_root.resolve(), a.campaign,
                     a.destination.resolve(), a.inventory.resolve())
    print(json.dumps(result))
