"""Stage provenance-only electronic parents for subsequent archived workers.

The original plan and independent review remain unchanged. Only source-map path
keys in a derived plan are relocated, with both original and derived bytes bound
explicitly. No electronic value, chemical field or numerical setting changes.
"""
import argparse
import json
from pathlib import Path
import shutil

from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json
from redox.validation.namespace_sources import NamespaceSources, read_checked_review


def relocate_plan(plan, reader, original_path):
    before = json.loads(original_path.read_text())
    if before != plan:
        raise ValueError('Only the exact original electronic plan may be relocated')
    result = dict(plan)
    for field in ('producer_sha256', 'input_sha256'):
        reader.verify_all(plan[field])
        result[field] = {str(reader.path(path)): digest for path, digest in plan[field].items()}
    result['provenance_path_relocation'] = dict(original_plan=str(original_path),
        original_plan_sha256=sha256_file(original_path), numerical_or_chemical_fields_changed=False,
        permitted_changes=['producer_sha256 path keys', 'input_sha256 path keys', 'this provenance record'])
    if ({k: v for k, v in result.items() if k not in ('producer_sha256', 'input_sha256', 'provenance_path_relocation')}
            != {k: v for k, v in plan.items() if k not in ('producer_sha256', 'input_sha256')}):
        raise ValueError('Plan relocation changed a chemical or numerical field')
    reader.recheck()
    return result


def stage(namespace, destination):
    report, binding = read_checked_review(namespace, 'energy')
    spec = json.loads(namespace.read_text())
    if (destination.exists() or destination.resolve().is_relative_to(Path(spec['prior_root']).resolve())
            or destination.resolve().is_relative_to(Path(spec['archive']).resolve())):
        raise ValueError('A fresh parent staging location outside source checkouts is required')
    case = Path(spec['attempt'])
    plan_path = case/'plan.json'
    original_hash = sha256_file(plan_path)
    reader = NamespaceSources(spec)
    plan = json.loads(plan_path.read_text())
    relocated = relocate_plan(plan, reader, plan_path)
    reader.verify_all(report['source_sha256'])
    destination.mkdir(parents=True)
    shutil.copyfile(plan_path, destination/'original_plan.json')
    shutil.copyfile(case/'input.xyz', destination/'input.xyz')
    atomic_write_json(destination/'plan.json', relocated)
    mapped = {str(reader.path(path)): digest for path, digest in report['source_sha256'].items()}
    mapped.update({str(Path(spec['original_campaign'])/'plan.json'): sha256_file(destination/'plan.json'),
                   str(Path(spec['original_raw_root'])/'experimental_energy.json'):
                       sha256_file(reader.path(Path(spec['original_raw_root'])/'experimental_energy.json')),
                   str(destination/'input.xyz'): sha256_file(destination/'input.xyz'),
                   str(destination/'original_plan.json'): original_hash,
                   binding['review']: binding['review_sha256']})
    derived = dict(report, source_sha256=mapped,
                   provenance_only_derivative_of_checked_review=True,
                   independent_recalculation_performed=False,
                   original_checked_review=binding)
    atomic_write_json(destination/'independent_energy_review.json', derived)
    staged_spec = dict(spec, attempt=str(destination), plan_sha256=sha256_file(destination/'plan.json'),
                       status='checked_energy_parent_with_explicit_source_path_relocation',
                       original_energy_namespace=str(namespace), original_energy_namespace_sha256=sha256_file(namespace))
    staged_path = destination/'namespace.json'
    atomic_write_json(staged_path, staged_spec)
    checked, staged_binding = read_checked_review(staged_path, 'energy')
    if any(checked[name] != report[name] for name in ('task', 'E_smd_eV', 'E_gas_eV', 'experimental_model', 'portable_method_id')):
        raise ValueError('Staged parent changed electronic identity or values')
    reader.recheck()
    if sha256_file(plan_path) != original_hash:
        raise ValueError('Original electronic plan changed during staging')
    atomic_write_json(destination/'derivation.json', dict(status='provenance_only_energy_parent_checked',
                     original_binding=binding, staged_binding=staged_binding,
                     original_plan_sha256=original_hash, derived_plan_sha256=staged_spec['plan_sha256'],
                     producer_sha256=sha256_file(Path(__file__)), numerical_or_chemical_fields_changed=False))
    return staged_spec


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--namespace', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    args = p.parse_args()
    from redox.workflows.resources import allocation_guard, cpu_allocation, host_capacity
    import os
    with allocation_guard():
        host_capacity(1, minimum_ram_gib=4)
        os.sched_setaffinity(0, cpu_allocation(1))
    print(json.dumps(stage(args.namespace.resolve(), args.out.resolve())))
