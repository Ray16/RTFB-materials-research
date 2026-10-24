#!/usr/bin/env python
"""Audit a fresh prefix and rerun existing singlet/doublet quantum controls.

Only dependency locations are rebound, to byte-identical qualified xTB artifacts.
The original optimizer, active energy, native Hessian and partition gates run in
the exact archived source checkout. No production record is replaced.
"""
import argparse
import hashlib
import importlib.metadata as metadata
import importlib.util
import json
import os
from pathlib import Path
import sys
import time


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write(path, value):
    temporary = path.with_name(path.name + '.tmp')
    temporary.write_text(json.dumps(value, indent=2) + '\n')
    temporary.replace(path)


def environment(restoration):
    report = json.loads((restoration / 'status.json').read_text())
    if report['status'] != 'fresh_environment_installed_GPU_qualification_pending':
        raise ValueError('Fresh installation must finish successfully before qualification')
    prefix, source = Path(report['prefix']).resolve(), Path(report['source_checkout']).resolve()
    if Path(sys.prefix).resolve() != prefix or not sys.flags.no_user_site:
        raise ValueError('Use the fresh interpreter with user site disabled')
    from redox.core.common import ROOT
    if ROOT.resolve() != source:
        raise ValueError('Project import escaped the exact archived source checkout')
    runtime = json.loads((source / 'environments/quantum_runtime.json').read_text())
    selected_lock = Path(report.get('selected_pip_lock', source / 'environments/quantum-pip.lock.txt'))
    if sha(selected_lock) != report['pip_lock_sha256']:
        raise ValueError('Selected installation lock changed')
    all_versions = {}
    for line in selected_lock.read_text().splitlines():
        if not line.strip() or line.lstrip().startswith('#'):
            continue
        name, expected = line.split('==', 1)
        actual = metadata.version(name)
        if actual != expected:
            raise ValueError('A selected distribution differs from its installation lock: ' + name)
        all_versions[name] = actual
    versions = {}
    for name, expected in runtime['versions'].items():
        actual = metadata.version(name)
        if actual != expected:
            raise ValueError('Recorded quantum version changed: ' + name)
        distribution = metadata.distribution(name)
        location = Path(distribution.locate_file('')).resolve()
        if not location.is_relative_to(prefix):
            raise ValueError('Scientific distribution escaped fresh prefix: ' + name)
        versions[name] = dict(version=actual, location=str(location))
    origins = {}
    for name in ('rdkit', 'numpy', 'scipy', 'ase', 'geometric', 'pyscf', 'gpu4pyscf', 'torch', 'cupy'):
        spec = importlib.util.find_spec(name)
        if spec is None or not Path(spec.origin).resolve().is_relative_to(prefix):
            raise ValueError('Scientific import escaped fresh prefix: ' + name)
        origins[name] = spec.origin
    library, binary = prefix / 'lib/libxtb.so.6.7.1', prefix / 'bin/xtb'
    if sha(library) != runtime['libxtb_sha256']:
        raise ValueError('Restored xTB library is not byte-identical to the qualified build')
    return dict(checked=time.time(), prefix=str(prefix), source=str(source),
        source_commit=report['source_commit'], versions=versions, module_origins=origins,
        selected_pip_lock_sha256=sha(selected_lock), all_distribution_versions=all_versions,
        libxtb=str(library), libxtb_sha256=sha(library), xtb=str(binary), xtb_sha256=sha(binary),
        qualification_driver_sha256=sha(__file__), restoration_status_sha256=sha(restoration / 'status.json'),
        user_site_disabled=True, source_archive_sha256=report['source_archive_sha256'])


def control_module(source):
    path = source / 'scripts/validation/quantum_restart_equivalence.py'
    spec = importlib.util.spec_from_file_location('restored_standard_quantum_controls', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def bind_paths(proof):
    from redox.qm import qrrho, native_thermal, projected_thermal
    if proof['libxtb_sha256'] != qrrho.LIBRARY_SHA256:
        raise ValueError('Runtime qualification hash differs from installed artifact')
    qrrho.LIBRARY = Path(proof['libxtb'])
    qrrho._checked_library.cache_clear()
    qrrho._checked_library()
    native_thermal.XTB = projected_thermal.XTB = Path(proof['xtb'])


def prepare(restoration, campaign, reference):
    proof = environment(restoration)
    module = control_module(Path(proof['source']))
    reference_data = json.loads(reference.read_text())
    if reference_data['status'] != 'original_environment_baseline_controls_independently_exported':
        raise ValueError('Independently admitted original control snapshot required')
    module.api.verify_sources(reference_data['source_sha256'])
    module.prepare(campaign)
    plan = json.loads((campaign / 'plan.json').read_text())
    plan['tasks'] = [t for t in plan['tasks'] if t['phase'] == 'baseline']
    plan['raw_root'] = str(campaign / 'raw')
    plan['role'] = 'Fresh-environment reconstruction of the two standard quantum baseline controls; no new restart claim'
    plan['producer_sha256'][str(Path(__file__).resolve())] = sha(__file__)
    plan['restored_environment'] = proof
    plan['reference_control_snapshot'] = str(reference)
    plan['input_sha256'][str(reference)] = sha(reference)
    plan['input_sha256'].update(reference_data['source_sha256'])
    write(campaign / 'plan.json', plan)
    write(campaign / 'environment_audit.json', proof)
    return proof


def worker(restoration, campaign, case):
    from redox.workflows.resources import verify_gate_owner
    verify_gate_owner()
    proof = environment(restoration)
    plan = json.loads((campaign / 'plan.json').read_text())
    if proof['source_commit'] != plan['restored_environment']['source_commit']:
        raise ValueError('Restored source changed after preparation')
    # Explicit process-local dependency configuration, recorded below. Scientific
    # sources, certificates, thresholds and library bytes are unchanged.
    bind_paths(proof)
    from redox.qm import qrrho
    module = control_module(Path(proof['source']))
    module.run_task(campaign, case + '__baseline')
    kernel = qrrho.kernel_execution()
    if Path(kernel['interpreter']).resolve() != Path(sys.executable).resolve():
        raise ValueError('Thermodynamic helper escaped the fresh interpreter')
    mapped = Path('/proc/self/maps').read_text()
    original_prefix = '/nfs/lambda_stor_01/homes/rzhu/miniforge3/envs/redox/'
    if original_prefix in mapped:
        raise ValueError('A shared library was loaded from the original quantum environment')
    proof.update(status='fresh_quantum_control_and_normal_native_admission_passed',
        case=case, kernel_execution=kernel, library_binding='explicit byte-identical fresh-prefix xTB paths',
        original_environment_shared_library_loaded=False, finished=time.time())
    write(campaign / (case + '_deployment.json'), proof)
    return proof


def review(restoration, campaign, reference):
    proof = environment(restoration)
    bind_paths(proof)
    from ase.io import read
    from ase import Atoms
    from redox.core.protocol import load_record
    from redox.qm.native_thermal import load_native_thermal
    from redox.workflows import polymer_campaign as api
    module = control_module(Path(proof['source']))
    plan = api.read_plan(campaign)
    api.verify_sources(plan['producer_sha256'])
    api.verify_sources(plan['input_sha256'])
    if str(reference) != plan['reference_control_snapshot']:
        raise ValueError('Comparison reference differs from the predeclared snapshot')
    previous = json.loads(reference.read_text())
    if previous['status'] != 'original_environment_baseline_controls_independently_exported':
        raise ValueError('Original controls lack independent admission')
    api.verify_sources(previous['source_sha256'])
    old_records = {r['case']: r for r in previous['records']}
    rows, witnesses = [], {}
    for case, definition in module.CASES.items():
        original = old_records[case]
        if original['mult'] != definition['mult'] or original['status'] != 'normal_primary_electronic_and_native_admission_passed':
            raise ValueError('Original comparison control identity/admission differs')
        records = [(Atoms(original['symbols'], positions=original['coordinates_A']),
            original['energy'], original['thermal'], dict(G_eV=original['G_eV']))]
        for root in (campaign,):
            path = Path(api.read_plan(root)['raw_root']) / (case + '__baseline/status.json')
            status = json.loads(path.read_text())
            if status['status'] != 'quantum_control_verified':
                raise ValueError('Every reference and fresh control must pass normal physical QC')
            api.verify_sources(status['output_sha256'])
            sd = Path(status['state_dir'])
            energy = load_record(sd, 0, definition['mult'])
            thermal = load_native_thermal(sd, 0, definition['mult'])
            if energy is None or thermal is None:
                raise ValueError('Backing electronic/native record no longer qualifies')
            if abs(energy['e_smd_eV'] + thermal['g_thermal_total_eV'] - status['G_eV']) > 1e-9:
                raise ValueError('Control G does not reconstruct')
            records.append((read(sd / 'opt.xyz'), energy, thermal, status))
            witnesses[str(path)] = sha(path)
            witnesses.update(status['output_sha256'])
        differences = module.restart_differences(*records)
        cutoff_error = max(abs(records[0][2]['nuclear_by_cutoff_eV'][k]
            - records[1][2]['nuclear_by_cutoff_eV'][k]) for k in ('25', '50', '100'))
        if cutoff_error > module.TOLERANCES['free_energy_eV']:
            raise ValueError('Nuclear partition differs beyond the predeclared G tolerance')
        deployment = campaign / (case + '_deployment.json')
        observed = json.loads(deployment.read_text())
        if observed['status'] != 'fresh_quantum_control_and_normal_native_admission_passed':
            raise ValueError('Fresh-prefix deployment isolation was not checked')
        witnesses[str(deployment)] = sha(deployment)
        rows.append(dict(case=case, mult=definition['mult'], differences=differences,
            maximum_cutoff_nuclear_difference_eV=cutoff_error))
    result = dict(status='fresh_compatible_environment_two_spin_controls_verified',
        checked=time.time(), controls=rows, tolerances=module.TOLERANCES,
        environment=proof, backing_sha256=witnesses,
        limits=['Two rigid controls do not establish universal chemical accuracy or flexible-polymer convergence',
            'Explicit fresh-prefix dependency bindings retain the xTB hash; shared-cluster gate configuration remains external',
            'The compatible CUDA runtime profile is distinct from the source-frozen recorded environment'])
    write(campaign / 'summary.json', result)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('prepare', 'worker', 'review'))
    parser.add_argument('--restoration', type=Path, required=True)
    parser.add_argument('--campaign', type=Path, required=True)
    parser.add_argument('--case', choices=('formaldehyde', 'methyl_radical'))
    parser.add_argument('--reference', type=Path)
    args = parser.parse_args()
    if args.action == 'worker' and not args.case:
        parser.error('--case is required for worker')
    if args.action in ('prepare', 'review') and not args.reference:
        parser.error('--reference is required for prepare/review')
    if args.action == 'prepare':
        result = prepare(args.restoration.resolve(), args.campaign.resolve(), args.reference.resolve())
    elif args.action == 'review':
        result = review(args.restoration.resolve(), args.campaign.resolve(), args.reference.resolve())
    else:
        result = worker(args.restoration.resolve(), args.campaign.resolve(), args.case)
    print(json.dumps(result))


if __name__ == '__main__':
    main()
