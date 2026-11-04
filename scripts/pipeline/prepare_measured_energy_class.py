#!/usr/bin/env python3
"""Prepare same-constitution energy successors inside the frozen source view.

This CPU-only driver calls the original electronic preparer and its resource
compatibility checks. Each geometry keeps its own independent review and exact
coordinates; only the measured resource budget transfers across conformers.
"""
import argparse
from contextlib import contextmanager
import json
import os
from pathlib import Path
import shlex
import subprocess
import time

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api

STOPPED_OWNERS = {}


class OwnershipUnverified(RuntimeError):
    pass


def owner_alive(owner):
    """Retain exact PID/start checks without assuming a moved remote Python."""
    host = owner['host'].split('.')[0]
    identity = host, int(owner['pid']), owner['process_start']
    if identity in STOPPED_OWNERS:
        return False
    path = Path('/proc')/str(int(owner['pid']))/'stat'
    if host == os.uname().nodename.split('.')[0]:
        alive = path.exists() and path.read_text().rsplit(')', 1)[1].split()[19] == owner['process_start']
        if not alive:
            STOPPED_OWNERS[identity] = dict(owner=owner, stopped=True)
        return alive
    code = ('from pathlib import Path; p=Path('+repr(str(path))+'); '
            'print(p.exists() and p.read_text().rsplit(")",1)[1].split()[19]=='+repr(owner['process_start'])+')')
    try:
        result = subprocess.run(['ssh', '-n', '-S', 'none', '-o', 'BatchMode=yes', '-o', 'ConnectTimeout=8', host,
                                 shlex.join(['/usr/bin/python3', '-I', '-S', '-c', code])],
                                capture_output=True, text=True, timeout=20, check=True)
    except (OSError, subprocess.SubprocessError) as error:
        raise OwnershipUnverified('Cannot verify original geometry owner') from error
    if result.stdout.strip() not in ('True', 'False'):
        raise OwnershipUnverified('Unexpected original geometry ownership response')
    alive = result.stdout.strip() == 'True'
    if not alive:
        STOPPED_OWNERS[identity] = dict(owner=owner, stopped=True)
    return alive


@contextmanager
def ownership_transport():
    original = api.load_module

    def load(name, path):
        module = original(name, path)
        if path == ROOT/'scripts/pipeline/run_polymer_spacing.py':
            module.owner_alive = owner_alive
        return module

    api.load_module = load
    try:
        yield
    finally:
        api.load_module = original


def ready(helper, entry, manifest):
    try:
        with ownership_transport():
            return helper.prepare_ready(entry, manifest)
    except OwnershipUnverified:
        return 'upstream_ownership_unverified'


def prepare_recovered(cfg, configuration, helper, out):
    """Use independently checked recoveries with unchanged declared identities."""
    rows, keys = [], set()
    for index, source in enumerate(cfg['geometry_subset_inputs']):
        original_path = Path(source['batch'])/'batch.json'
        original = json.loads(original_path.read_text())
        api.verify_sources(original['source_sha256'])
        requested = set(source['keys'])
        selected = [entry for entry in original['entries'] if entry['key'] in requested]
        if not requested or len(selected) != len(requested):
            raise ValueError('Recovery subset must name exact declared keys')
        for entry in selected:
            before = api.read_plan(Path(cfg['recovered_from'][entry['key']]))
            after = api.read_plan(Path(entry['campaign']))
            if before['tasks'] != after['tasks'] or entry['key'] in keys:
                raise ValueError('Recovery changed the declared task or repeated its identity')
            keys.add(entry['key'])
        subset = out/('recovery_subset_'+str(index))
        subset.mkdir()
        atomic_write_json(subset/'batch.json', dict(original, entries=selected,
                          source_sha256={**original['source_sha256'], str(original_path): api.digest(original_path)}))
        target = out/(out.name+'__recovery_'+str(index))
        helper.declare(subset, Path(source['reviews']), Path(cfg['reference']), Path(cfg['resource_review']), target)
        manifest = json.loads((target/'batch.json').read_text())
        manifest['maximum_RSS_GiB'] = cfg['maximum_RSS_GiB']
        manifest['source_sha256'].update({str(configuration): api.digest(configuration),
                                         str(Path(__file__).resolve()): api.digest(Path(__file__))})
        for entry in manifest['entries']:
            entry['geometry_review'] = str(Path(source['reviews'])/(Path(source['batch']).name+'__'+entry['key']+'.json'))
        atomic_write_json(target/'batch.json', manifest)
        for entry in manifest['entries']:
            rows.append(dict(key=entry['key'], campaign=entry['campaign'], status=ready(helper, entry, manifest),
                             batch=str(target/'batch.json'), geometry_review=entry['geometry_review'],
                             recovered_from=cfg['recovered_from'][entry['key']]))
    return dict(status='same_constitution_energy_dependencies_prepared', entries=rows,
                declared_states=len(keys), maximum_RSS_GiB=cfg['maximum_RSS_GiB'],
                primary_candidate_admission=False, quantum_jobs_submitted=0, native_G_computed=False,
                stopped_geometry_owners=list(STOPPED_OWNERS.values()))


def prepare(configuration):
    cfg = json.loads(configuration.read_text())
    if Path(cfg['prior_root']) != ROOT or cfg['primary_candidate_admission'] is not False:
        raise ValueError('Explicit original source view and method-control scope required')
    api.verify_sources(cfg['source_sha256'])
    helper = api.load_module('measured_class_original_preparer', ROOT/'scripts/pipeline/prepare_smooth_electronic_batch.py')
    point = api.read_plan(Path(cfg['reference']))
    resource = json.loads(Path(cfg['resource_review']).read_text())
    helper.resource_compatible(point, point, resource)
    api.verify_sources(resource['source_sha256'])
    bound = cfg['maximum_RSS_GiB']
    if not 1.25*resource['peak_host_RSS_bytes']/2**30 <= bound <= point['maximum_host_RSS_GiB']:
        raise ValueError('Measured memory headroom and original pilot ceiling required')
    out = Path(cfg['logical_output'])
    out.mkdir(parents=True, exist_ok=False)
    if cfg.get('geometry_subset_inputs'):
        result = prepare_recovered(cfg, configuration, helper, out)
        atomic_write_json(out/'preparation.json', result)
        return result
    batches = list(cfg['geometry_batches'])
    pilot = Path(cfg['pilot_batch'])
    original = json.loads((pilot/'batch.json').read_text())
    api.verify_sources(original['source_sha256'])
    matching = []
    for entry in original['entries']:
        plan = api.read_plan(Path(entry['campaign']))
        try:
            helper.resource_compatible(plan, point, resource)
        except ValueError:
            continue
        matching.append(entry)
    if len(matching) != 2:
        raise ValueError('Exactly two explicit configuration pilots required')
    subset = out/'pilot_geometry_subset'
    subset.mkdir()
    atomic_write_json(subset/'batch.json', dict(original, entries=matching,
        source_sha256={**original['source_sha256'], str(pilot/'batch.json'): api.digest(pilot/'batch.json')}))
    batches.append(dict(batch=str(subset), reviews=cfg['pilot_reviews'], original_review_prefix=pilot.name))
    statuses, identities = [], set()
    for index, source in enumerate(batches):
        batch, reviews = Path(source['batch']), Path(source['reviews'])
        target = out/(out.name+'__'+str(index)+'__'+batch.name)
        helper.declare(batch, reviews, Path(cfg['reference']), Path(cfg['resource_review']), target)
        path = target/'batch.json'
        manifest = json.loads(path.read_text())
        manifest.update(maximum_RSS_GiB=bound,
                        resource_budget_derivation=dict(measured_peak_GiB=resource['peak_host_RSS_bytes']/2**30,
                                                       guard_GiB=bound, numerical_protocol_changed=False))
        manifest['source_sha256'].update({str(configuration): api.digest(configuration),
                                         str(Path(__file__).resolve()): api.digest(Path(__file__))})
        for entry in manifest['entries']:
            if 'original_review_prefix' in source:
                entry['geometry_review'] = str(reviews/(source['original_review_prefix']+'__'+entry['key']+'.json'))
            override = cfg.get('geometry_review_overrides', {}).get(entry['key'])
            if override:
                checked = json.loads(Path(override).read_text())
                if checked['campaign'] != entry['geometry_campaign']:
                    raise ValueError('Geometry-review override belongs to a different exact input')
                entry['geometry_review'] = override
            if entry['key'] in identities:
                raise ValueError('Conformer identity repeated in declared geometry batches')
            identities.add(entry['key'])
        atomic_write_json(path, manifest)
        for entry in manifest['entries']:
            state = ready(helper, entry, manifest)
            statuses.append(dict(key=entry['key'], campaign=entry['campaign'], status=state,
                                 batch=str(path), geometry_review=entry['geometry_review']))
            atomic_write_json(out/'preparation_progress.json', dict(entries=statuses,
                              primary_candidate_admission=False, quantum_jobs_submitted=0))
    if len(identities) != 40:
        raise ValueError('The two declared twenty-input configurations must total forty states')
    api.verify_sources(cfg['source_sha256'])
    result = dict(status='same_constitution_energy_dependencies_prepared', entries=statuses,
                  declared_states=40, maximum_RSS_GiB=bound, primary_candidate_admission=False,
                  quantum_jobs_submitted=0, native_G_computed=False,
                  stopped_geometry_owners=list(STOPPED_OWNERS.values()))
    atomic_write_json(out/'preparation.json', result)
    return result


def refresh(configuration):
    """Recheck held dependencies without replacing any prepared or live state."""
    cfg = json.loads(configuration.read_text())
    api.verify_sources(cfg['source_sha256'])
    out = Path(cfg['logical_output'])
    previous = json.loads((out/'preparation.json').read_text())
    history = out/'preparation_history'
    history.mkdir(exist_ok=True)
    (history/(str(time.time_ns())+'.json')).write_bytes((out/'preparation.json').read_bytes())
    helper = api.load_module('refresh_original_energy_dependency_preparer', ROOT/'scripts/pipeline/prepare_smooth_electronic_batch.py')
    rows = []
    for batch_path in dict.fromkeys(row['batch'] for row in previous['entries']):
        manifest = json.loads(Path(batch_path).read_text())
        api.verify_sources(manifest['source_sha256'])
        for entry in manifest['entries']:
            state = ready(helper, entry, manifest)
            rows.append(dict(key=entry['key'], campaign=entry['campaign'], status=state,
                             batch=batch_path, geometry_review=entry['geometry_review']))
    result = dict(previous, entries=rows, stopped_geometry_owners=list(STOPPED_OWNERS.values()),
                  dependency_refresh_source=str(Path(__file__).resolve()),
                  dependency_refresh_source_sha256=api.digest(Path(__file__)))
    atomic_write_json(out/'preparation.json', result)
    return result


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--configuration', type=Path, required=True)
    p.add_argument('--refresh', action='store_true')
    args = p.parse_args()
    result = (refresh if args.refresh else prepare)(args.configuration.resolve())
    print(json.dumps(dict(status=result['status'], declared_states=result['declared_states'])))
