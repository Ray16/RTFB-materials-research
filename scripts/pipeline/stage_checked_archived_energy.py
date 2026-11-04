#!/usr/bin/env python3
"""Stage byte-identical reviewed electronic evidence for a native handoff."""
import argparse
import json
from pathlib import Path, PurePosixPath
import shutil

from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json
from redox.validation.archived_sources import ArchivedSources
from redox.validation.namespace_sources import read_checked_review


def stage(archive, prior, review_path, destination):
    report = json.loads(review_path.read_text())
    if (report['status'] != 'alternate_electronic_geometry_and_resource_checked_native_G_pending'
            or report['primary_candidate_admission'] is not False
            or report['native_Hessian_or_G_qualified'] is not False):
        raise ValueError('A completed independent electronic review without native G is required')
    reader = ArchivedSources(prior, archive, relative_to_prior=True)
    reader.verify_all(report['source_sha256'])
    original_case = Path(report['campaign'])
    plan = json.loads(reader.path(original_case/'plan.json').read_text())
    if plan['tasks'] != [report['task']] or plan['experimental_model'] != report['experimental_model']:
        raise ValueError('Exact electronic task/model differs')
    if destination.exists() or destination.is_relative_to(archive) or destination.is_relative_to(prior):
        raise ValueError('Use a fresh staged-evidence location outside source checkouts')
    case, raw = destination/'campaign', destination/'raw'
    case.mkdir(parents=True)
    raw.mkdir()
    original_raw = Path(plan['raw_root'])
    original_raw_parent = original_raw.parent
    copies = {}
    for path, expected in report['source_sha256'].items():
        original = Path(path) if Path(path).is_absolute() else prior/path
        target = None
        for prefix, replacement in ((original_case, case), (original_raw_parent, raw)):
            try:
                target = replacement/original.relative_to(prefix)
                break
            except ValueError:
                continue
        if target is None:
            continue
        target.parent.mkdir(parents=True, exist_ok=True)
        source = reader.verify(path, expected)
        shutil.copyfile(source, target)
        if sha256_file(target) != expected:
            raise ValueError('Staging changed a reviewed input')
        copies[str(target)] = expected
    for name in ('plan.json', 'input.xyz'):
        source = reader.path(original_case/name)
        target = case/name
        if target.exists() and sha256_file(target) != sha256_file(source):
            raise ValueError('Conflicting staged electronic input')
        if not target.exists():
            shutil.copyfile(source, target)
        copies[str(target)] = sha256_file(target)
    shutil.copyfile(review_path, case/'independent_energy_review.json')
    spec = dict(status='byte_identical_archived_electronic_evidence_staged', archive=str(archive), prior_root=str(prior),
                original_campaign=str(original_case), attempt=str(case), original_raw_parent=str(original_raw_parent),
                original_raw_root=str(original_raw), raw_parent=str(raw), key=report['task']['key'],
                method_id=report['portable_method_id'], plan_sha256=sha256_file(case/'plan.json'),
                geometry_sha256=plan['scored_geometry_sha256'], original_independent_review=str(review_path),
                original_independent_review_sha256=sha256_file(review_path),
                quantum_calculations_run=False, numerical_or_chemical_fields_changed=False)
    path = destination/'namespace.json'
    atomic_write_json(path, spec)
    _, binding = read_checked_review(path, 'energy')
    reader.recheck()
    atomic_write_json(destination/'staging.json', dict(status='staged_electronic_evidence_bindings_checked',
                     copied_files_sha256=copies, checked_binding=binding,
                     original_review_sha256=sha256_file(review_path), quantum_calculations_run=False))
    return spec


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('archive', 'prior-root', 'review', 'out'):
        p.add_argument('--'+name, type=Path, required=True)
    a = p.parse_args()
    print(json.dumps(stage(a.archive.resolve(), a.prior_root.resolve(), a.review.resolve(), a.out.resolve())))
