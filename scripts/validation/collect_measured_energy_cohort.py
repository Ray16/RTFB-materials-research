#!/usr/bin/env python3
"""Audit a complete declared energy roster with explicit recovery precedence."""
import argparse
from collections import Counter, defaultdict, deque
import json
from pathlib import Path
import re
import shutil

from redox.core.common import ROOT
from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json
from redox.validation.namespace_sources import read_checked_review
from redox.workflows.archived_native_namespace import verify_spec


def collect(manifests, output):
    if output.exists():
        raise FileExistsError('Preserve previous cohort snapshots')
    selected, declared, producers = {}, set(), {str(Path(__file__).resolve()): sha256_file(__file__)}
    copied = {}
    for path in manifests:
        data = json.loads(path.read_text())
        if data['primary_candidate_admission'] or data['native_G_computed'] or data['quantum_jobs_submitted']:
            raise ValueError('Only exact prepared method-control dependencies supported')
        producers[str(path)] = sha256_file(path)
        copied[path.parent.name+'.json'] = path
        config = Path(data['source_configuration'])
        if sha256_file(config) != data['source_configuration_sha256']:
            raise ValueError('Preparation configuration changed')
        producers[str(config)] = sha256_file(config)
        namespace_by_key = {}
        for name in data['namespaces']:
            spec = verify_spec(Path(name))
            if spec['key'] in namespace_by_key:
                raise ValueError('Prepared namespace repeats a conformer identity')
            namespace_by_key[spec['key']] = (name, spec)
        for entry in data['entries']:
            key = entry['key']
            declared.add(key)
            if entry['status'] not in ('prepared', 'existing_checked_reference'):
                continue
            if key in selected:
                raise ValueError('Two eligible versions of one declared conformer; resolve explicitly')
            if entry['status'] == 'prepared':
                name, spec = namespace_by_key[key]
                plan = json.loads((Path(spec['attempt'])/'plan.json').read_text())
                row = dict(key=key, status='electronic_prepared_not_computed', namespace=name,
                           namespace_sha256=sha256_file(name), task=plan['tasks'][0],
                           geometry_sha256=plan['scored_geometry_sha256'], plan_sha256=spec['plan_sha256'],
                           maximum_host_RSS_GiB=plan['maximum_host_RSS_GiB'],
                           method_id=plan['experimental_method_id'], source_geometry=plan['source_geometry_campaign'])
            else:
                cfg = json.loads(config.read_text())
                report, binding = read_checked_review(Path(cfg['energy_namespace']), 'energy')
                if report['task']['key'] != key or report['campaign'] != json.loads(Path(cfg['energy_namespace']).read_text())['original_campaign']:
                    raise ValueError('Cached resource point has another exact identity')
                row = dict(key=key, status='electronic_independently_checked_reference', task=report['task'],
                           method_id=report['portable_method_id'], binding=binding,
                           E_smd_eV=report['E_smd_eV'], E_gas_eV=report['E_gas_eV'])
            if row['task']['key'] != key or row['task']['candidate_id'] != 'exp_phthalimide_control':
                raise ValueError('The declared graph or method-control identity differs')
            row.update(source_manifest=str(path), source_manifest_sha256=producers[str(path)],
                       recovery_origin=entry.get('recovered_from'), native_G_computed=False,
                       primary_candidate_admission=False)
            selected[key] = row
    if set(selected) != declared:
        raise ValueError('Declared inputs still lack exact prepared or checked replacements: '+str(sorted(declared-set(selected))))
    counts = Counter((r['task']['stereoisomer'], tuple(r['task']['occupied_sites'])) for r in selected.values())
    expected = {(configuration, sites) for configuration in ('all_same', 'alternating')
                for sites in ((3,), (4,), (2, 3), (2, 4))}
    if len(selected) != 160 or set(counts) != expected or set(counts.values()) != {20}:
        raise ValueError('The four remaining occupancy classes need two twenty-input configurations each')
    for signature in expected:
        group = [r['task'] for r in selected.values()
                 if (r['task']['stereoisomer'], tuple(r['task']['occupied_sites'])) == signature]
        if (len({r['graph_sha256'] for r in group}) != 1
                or len({json.dumps(r['backbone_stereo_fingerprint']) for r in group}) != 1):
            raise ValueError('Graph or backbone stereochemistry differs within a sampling group')
    if len({r['method_id'] for r in selected.values()}) != 1:
        raise ValueError('Different solvent methods cannot share this cohort')
    output.mkdir(parents=True)
    (output/'inputs').mkdir()
    for name, path in copied.items():
        shutil.copyfile(path, output/'inputs'/name)
    policy = ROOT.parent/'.config/gpu_reserve/config.env'
    hosts = re.search(r'^ALLOWED_HOSTS="([^"]+)"', policy.read_text(), re.M)[1].split()
    groups = defaultdict(deque)
    for key in sorted(selected):
        row = selected[key]
        if row['status'] == 'electronic_prepared_not_computed':
            groups[tuple(row['task']['occupied_sites'])].append(row['namespace'])
    order = []
    while any(groups.values()):
        for key in ((2, 4), (2, 3), (4,), (3,)):
            if groups[key]:
                order.append(groups[key].popleft())
    queue = dict(specs=order, hosts=hosts, declared_upstream_states=160,
                 checked_references=sum(r['status']=='electronic_independently_checked_reference' for r in selected.values()),
                 native_G_computed=False, source_policy_sha256=sha256_file(policy))
    atomic_write_json(output/'queue_inputs.json', queue)
    result = dict(status='declared_energy_dependencies_complete_quantum_evaluation_pending',
                  declared_states=160, prepared_states=len(order), checked_references=queue['checked_references'],
                  graph_configuration_counts=[dict(configuration=k[0], occupied_sites=list(k[1]), states=v)
                                              for k, v in sorted(counts.items())], records=list(selected.values()),
                  source_sha256=producers, native_G_computed=False, global_ensemble_convergence_proven=False,
                  primary_candidate_admission=False, scientific_accuracy_inferred_from_preparation=False)
    if any(sha256_file(path) != digest for path, digest in producers.items()):
        raise ValueError('Cohort inputs changed during capture')
    atomic_write_json(output/'summary.json', result)
    return result


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--manifest', type=Path, action='append', required=True)
    p.add_argument('--out', type=Path, required=True)
    a = p.parse_args()
    from redox.workflows.resources import allocation_guard, cpu_allocation, host_capacity
    import os
    with allocation_guard():
        host_capacity(1, minimum_ram_gib=4)
        os.sched_setaffinity(0, cpu_allocation(1))
    result = collect([path.resolve() for path in a.manifest], a.out.resolve())
    print(json.dumps({k: result[k] for k in ('status', 'declared_states', 'prepared_states', 'checked_references')}))
