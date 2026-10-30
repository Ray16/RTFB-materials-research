#!/usr/bin/env python3
"""Run a predeclared measured-class energy successor without relabeling a pilot.

This process-local policy adapter reuses the established namespace launcher. It
does not edit its files while other collectors use them or change QM equations.
"""
import argparse
from contextlib import contextmanager
import json
import math
from pathlib import Path
import shutil
import threading

from rdkit import Chem

from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json
from redox.validation.archived_sources import ArchivedSources
from redox.workflows import archived_energy_namespace as base


def validate_measured(plan, prior, archive, *, source_reader=None):
    bound = plan['maximum_host_RSS_GiB']
    if (len(plan['tasks']) != 1 or bound not in (160, 240)
            or plan.get('prospective_resource_pilot') is not False
            or plan['execution_driver'] != str(prior/'scripts/pipeline/smooth_polymer_energy_bounded.py')
            or plan['native_G_computed'] is not False or plan['primary_candidate_admission'] is not False
            or not plan['energy_only_execution'] or plan['maximum_peak_GPU_MiB'] != 31744
            or plan['experimental_model']['solvent'] != 'acetonitrile'):
        raise ValueError('Only the declared 160/240-GiB measured neutral-control classes are supported')
    task = plan['tasks'][0]
    if ((task['charge'], task['mult']) != (0, 1)
            or task.get('candidate_id') != 'exp_phthalimide_control'
            or task.get('role') != 'interaction_method_control; not synthesis admission'):
        raise ValueError('Exact neutral-singlet method-control scope required')
    reader = source_reader or ArchivedSources(prior, archive, relative_to_prior=True)
    if str(reader.prior) != str(prior) or reader.archive != archive.resolve():
        raise ValueError('Resource source reader has different declared checkout roots')
    reader.verify_all(plan['producer_sha256'])
    reader.verify_all(plan['input_sha256'])
    reference_path = plan['electronic_RAM_class_review']
    if reference_path not in plan['input_sha256']:
        raise ValueError('Measured resource review is not bound by the exact plan')
    witness = json.loads(reader.verify(reference_path, plan['input_sha256'][reference_path]).read_text())
    reader.verify_all(witness['source_sha256'])
    point_path = str(Path(plan['measured_alternate_model_electronic_reference'])/'plan.json')
    if point_path not in plan['input_sha256']:
        raise ValueError('Measured reference plan is not bound by the successor')
    point = json.loads(reader.verify(point_path, plan['input_sha256'][point_path]).read_text())
    if (any(plan[name] != point[name] for name in ('quantum_settings', 'runtime_profile', 'experimental_model'))
            or not bound <= point['maximum_host_RSS_GiB']):
        raise ValueError('Reference numerical settings/runtime or resource ceiling differs')
    if (witness['status'] != 'alternate_electronic_geometry_and_resource_checked_native_G_pending'
            or witness['native_Hessian_or_G_qualified'] or witness['primary_candidate_admission']
            or witness['portable_method_id'] != plan['experimental_method_id']
            or witness['experimental_model'] != plan['experimental_model']
            or witness['campaign'] != plan['measured_alternate_model_electronic_reference']
            or not 0 < witness['peak_host_RSS_bytes']*1.25 <= bound*2**30
            or not 0 < witness['observed_peak_GPU_MiB'] <= 31744):
        raise ValueError('Measured electronic resource witness differs or lacks required headroom')
    other = witness['task']
    if any(task[k] != other[k] for k in ('atoms', 'active_basis_functions', 'charge', 'mult')):
        raise ValueError('Resource witness has different state, atom count or basis size')
    molecules = [Chem.MolFromSmiles(t['smiles']) for t in (task, other)]
    if any(m is None for m in molecules) or len({Chem.MolToSmiles(m, isomericSmiles=False) for m in molecules}) != 1:
        raise ValueError('Resource witness has a different constitution')
    reader.recheck()
    return dict(resource_review=reference_path, resource_review_sha256=plan['input_sha256'][reference_path],
                measured_peak_GiB=witness['peak_host_RSS_bytes']/2**30, bound_GiB=bound,
                source_bindings=len(reader.receipts), resource_budget_only_no_energy_transferred=True)


@contextmanager
def measured_policy(archive):
    if threading.current_thread() is not threading.main_thread() or threading.active_count() != 1:
        raise RuntimeError('Install the isolated launcher policy before starting Python threads')
    old = base.validate_plan
    base.validate_plan = lambda plan, prior: validate_measured(plan, prior, archive)
    try:
        yield
    finally:
        base.validate_plan = old


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('action', choices=('prepare', 'probe', 'launch'))
    for name in ('archive', 'prior-root', 'campaign', 'destination', 'inventory', 'spec', 'out'):
        p.add_argument('--'+name, type=Path)
    p.add_argument('--gpu', type=int)
    a = p.parse_args()
    if a.action == 'prepare':
        if any(getattr(a, name) is None for name in ('archive', 'prior_root', 'campaign', 'destination', 'inventory')):
            p.error('Explicit source, exact campaign, stopped-owner inventory and fresh destination required')
        with measured_policy(a.archive.resolve()):
            spec = base.prepare(a.archive.resolve(), a.prior_root.resolve(), a.campaign,
                                a.destination.resolve(), a.inventory.resolve())
        spec.update(status='prepared_measured_class_energy_followup', prospective_RAM_experiment=False,
                    measured_resource_class_revalidated=True)
        atomic_write_json(a.destination.resolve()/'namespace.json', spec)
        print(json.dumps(spec))
    elif a.action == 'probe':
        print(json.dumps(base.probe(a.spec.resolve())))
    else:
        if a.spec is None or a.out is None or a.gpu is None:
            p.error('Explicit namespace, fresh placement directory and GPU required')
        spec = json.loads(a.spec.read_text())
        a.out.resolve().parent.mkdir(parents=True, exist_ok=True)
        snapshot = a.out.resolve().parent/'measured_validation_source.py'
        if snapshot.exists():
            raise FileExistsError('Preserve previous validation source snapshots')
        shutil.copyfile(Path(__file__), snapshot)
        with measured_policy(Path(spec['archive'])):
            child, record = base.launch(a.spec.resolve(), a.gpu, a.out.resolve())
        atomic_write_json(a.out.resolve()/'measured_policy.json',
                          dict(source=str(snapshot), source_sha256=sha256_file(snapshot),
                               energy_launcher_sha256=sha256_file(Path(base.__file__)),
                               scientific_settings_changed=False))
        print(json.dumps(record), flush=True)
        record.update(returncode=child.wait(), status='measured_energy_successor_returned_review_required')
        atomic_write_json(a.out.resolve()/'completion.json', record)


if __name__ == '__main__':
    main()
