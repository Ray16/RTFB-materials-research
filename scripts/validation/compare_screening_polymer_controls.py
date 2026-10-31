#!/usr/bin/env python3
"""Compare repeated installed/corrected polymer forces without causal promotion."""
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path

import numpy as np

from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json
from redox.workflows.molecular_native_namespace import verify


def compare(reviews):
    if len(reviews) != 4:
        raise ValueError('Two complete repeats per screening variant required')
    records = {}
    sources = {}
    namespaces = set()
    source_namespace = None
    point_sources = None
    for path in reviews:
        review = json.loads(path.read_text())
        if (review['status'] != 'three_SCF_gradient_diagnostic_reconstructed_no_native_admission'
                or not review['full_execution_qualified'] or review['full_native_G_qualified']):
            raise ValueError('A closed independently reviewed gradient control is required')
        spec = verify(Path(review['diagnostic_namespace']))
        if review['diagnostic_namespace'] in namespaces:
            raise ValueError('Repeating one namespace is not an independent repeat')
        namespaces.add(review['diagnostic_namespace'])
        case = Path(spec['case'])
        plan = json.loads((case / 'plan.json').read_text())
        receipt = case / 'screening_variant.json'
        variant = json.loads(receipt.read_text())['variant']
        if variant != plan['screening_variant'] or variant not in ('installed', 'corrected'):
            raise ValueError('Control implementation identity differs')
        if source_namespace is None:
            source_namespace = plan['original_namespace']
            point_sources = plan['points']
        if plan['original_namespace'] != source_namespace or plan['points'] != point_sources:
            raise ValueError('Different source geometries entered a repeated comparison')
        for name, expected in review['source_sha256'].items():
            if sha256_file(name) != expected:
                raise ValueError('Reviewed control backing changed')
        sources.update(review['source_sha256'])
        sources[str(path)] = sha256_file(path)
        sources[str(receipt)] = sha256_file(receipt)
        if variant == 'corrected':
            directory = case / 'screening_implementation'
            metadata = json.loads((directory / 'implementation.json').read_text())
            helper_path = Path(spec['sources']) / 'helper.py'
            definition = importlib.util.spec_from_file_location('review_exact_barrier_source', helper_path)
            helper = importlib.util.module_from_spec(definition)
            definition.loader.exec_module(helper)
            expected = helper.standalone_source((Path(spec['sources']) / 'kernel.cu').read_text(), metadata['constants'])
            if (metadata['status'] != 'run_local_screening_context_closed_qualification_separate'
                    or metadata['installed_library_files_modified'] or metadata['physics_model_changed']
                    or metadata['kernel_sha256'] != hashlib.sha256(expected.encode()).hexdigest()
                    or (directory / 'screening_kernel.cu').read_text() != expected):
                raise ValueError('Declared upstream correction source differs')
            for name in ('implementation.json', 'screening_kernel.cu', 'corrected_host.py', 'original_host.py'):
                sources[str(directory / name)] = sha256_file(directory / name)
        points = []
        for i in range(3):
            folder = case / f'point{i}'
            row = json.loads((folder / 'result.json').read_text())
            arrays = {name: np.load(folder / (name + '.npy'), allow_pickle=False)
                      for name in ('gradient', 'solute', 'solvent', 'cds')}
            points.append(dict(energy_Ha=row['E_Ha'], **arrays))
        records.setdefault(variant, []).append(dict(points=points,
            source_Hessian_column_relative_error=review['repeated_relative_error_against_source_Hessian']))
    if set(records) != {'installed', 'corrected'} or any(len(v) != 2 for v in records.values()):
        raise ValueError('Exactly two independent controls per variant required')
    differences = []
    for i in range(3):
        values = {variant: [r['points'][i] for r in rows] for variant, rows in records.items()}
        row = dict(source_evaluation=point_sources[i]['evaluation'], point=i, components={})
        for name in ('energy_Ha', 'gradient', 'solute', 'solvent', 'cds'):
            row['components'][name] = dict(
                installed_repeat_maximum=float(np.max(abs(np.asarray(values['installed'][1][name]) - values['installed'][0][name]))),
                corrected_repeat_maximum=float(np.max(abs(np.asarray(values['corrected'][1][name]) - values['corrected'][0][name]))),
                variant_maximum=float(max(np.max(abs(np.asarray(a[name]) - b[name]))
                    for a in values['installed'] for b in values['corrected'])))
        differences.append(row)
    if any(sha256_file(p) != h for p, h in sources.items()):
        raise ValueError('Source changed during comparison')
    return dict(status='two_repeats_per_screening_variant_polymer_forces_compared',
        original_namespace=source_namespace, points=differences,
        source_Hessian_column_relative_errors={v: [r['source_Hessian_column_relative_error'] for r in rows] for v, rows in records.items()},
        largest_pairwise_gradient_difference_Ha_Bohr=max(
            value for row in differences for value in row['components']['gradient'].values()),
        unique_cause_identified=False,
        limitation='Three source-fixed geometries and two repeats per variant only; no bound on all candidate or polymer states.',
        corrected_polymer_Hessian_qualified=False, native_G_qualified=False,
        production_implementation_adopted=False, source_sha256=sources, producer_sha256=sha256_file(__file__))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--review', type=Path, action='append', required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    from redox.workflows.resources import allocation_guard, host_capacity, cpu_allocation
    with allocation_guard():
        host_capacity(1, minimum_ram_gib=4)
        os.sched_setaffinity(0, cpu_allocation(1))
    if args.out.exists():
        raise FileExistsError('Preserve previous implementation comparisons')
    result = compare([p.resolve() for p in args.review])
    atomic_write_json(args.out, result)
    print(json.dumps({k: result[k] for k in ('status', 'source_Hessian_column_relative_errors', 'unique_cause_identified')}))
