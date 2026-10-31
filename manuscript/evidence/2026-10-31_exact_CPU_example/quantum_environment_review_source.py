#!/usr/bin/env python
"""Reconstruct fresh-deployment state sums and Hessian partitions for publication."""
import argparse
import json
import math
import os
from pathlib import Path
import time

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.core.protocol import load_record
from redox.workflows.polymer_campaign import digest, verify_sources
from redox.workflows.resources import allocation_guard, cpu_allocation, host_capacity


def publish(restoration, campaign, evidence, out):
    with allocation_guard():
        host_capacity(1)
        os.sched_setaffinity(0, cpu_allocation(1))
    import numpy as np
    from ase.io import read
    from redox.qm import qrrho
    from redox.qm.native_thermal import cache_path, NATIVE_THERMAL
    from redox.qm.projected_thermal import cache_path as partition_path
    report = json.loads((campaign / 'summary.json').read_text())
    if report['status'] != 'fresh_compatible_environment_two_spin_controls_verified':
        raise ValueError('Completed two-spin deployment comparison required')
    verify_sources(report['backing_sha256'])
    installation = json.loads((restoration / 'status.json').read_text())
    if not installation['dependency_check_passed']:
        raise ValueError('Installation dependency check did not pass')
    reference = json.loads((evidence / 'reference_control_snapshot.json').read_text())
    verify_sources(reference['source_sha256'])
    previous = {r['case']: r for r in reference['records']}
    plan = json.loads((campaign / 'plan.json').read_text())
    verify_sources(plan['producer_sha256'])
    verify_sources(plan['input_sha256'])
    states, comparisons = [], []
    for comparison in report['controls']:
        case, mult = comparison['case'], comparison['mult']
        status_path = Path(plan['raw_root']) / (case + '__baseline/status.json')
        status = json.loads(status_path.read_text())
        if status['status'] != 'quantum_control_verified':
            raise ValueError('Normal full quantum/native admission required')
        guard = evidence / 'complete/guards' / case
        owner = json.loads((guard / 'owner.json').read_text())
        residency = json.loads((guard / 'state.json').read_text())
        if (not owner.get('guard_stopped') or owner['pid'] != status['owner']['pid']
                or owner['process_start'] != status['owner']['process_start']
                or residency['foreign_resident_PIDs'] or not residency['own_residency_observed']
                or (guard / 'interruption.json').exists() or (guard / 'observation_error.json').exists()):
            raise ValueError('Complete exclusive-device ownership proof required')
        sd = Path(status['state_dir'])
        energy = load_record(sd, 0, mult)
        if energy is None:
            raise ValueError('Fresh active electronic parent failed independent admission')
        native_path, identity = cache_path(sd, 0, mult)
        native = json.loads(native_path.read_text())
        if (native['identity'] != identity or native['identity']['protocol'] != NATIVE_THERMAL
                or native['status'] != 'derivative_checked_solution_thermal'
                or not native['stationary_minimum']):
            raise ValueError('Native physical/identity record differs')
        target, ident = partition_path(native_path, sd / 'opt.xyz')
        thermal = json.loads(target.read_text())
        if thermal['identity'] != ident or thermal['status'] != 'projected_internal_qrrho_checked':
            raise ValueError('Projected partition is not bound to the stored Hessian/geometry')
        atoms = read(sd / 'opt.xyz')
        matrix = np.load(native_path.with_suffix('') / 'hessian_cartesian.npy', allow_pickle=False)
        rebuilt = qrrho.nuclear_thermals(matrix, atoms.positions, atoms.numbers, thermal['symmetry_number'])
        frequency_error = float(np.max(abs(np.asarray(rebuilt['frequencies_cm']) - thermal['frequencies_cm'])))
        nuclear_error = max(abs(rebuilt['nuclear_by_cutoff_eV'][k] - thermal['nuclear_by_cutoff_eV'][k])
            for k in ('25', '50', '100'))
        spin = -8.617333262145e-5 * 298.15 * math.log(mult)
        G = energy['e_smd_eV'] + rebuilt['nuclear_by_cutoff_eV']['50'] + spin
        if frequency_error > 1e-6 or nuclear_error > 1e-8 or abs(G - status['G_eV']) > 1e-9:
            raise ValueError('Fresh frequencies or free energy failed independent reconstruction')
        fresh = dict(case=case, mult=mult, symbols=atoms.get_chemical_symbols(),
            coordinates_A=atoms.positions.tolist(), energy={k:energy[k] for k in ('e_smd_eV', 'e_gas_eV')},
            thermal={k:thermal[k] for k in ('frequencies_cm', 'nuclear_by_cutoff_eV')},
            G_eV=G, status='normal_primary_electronic_and_native_admission_passed')
        states.append(dict(case=case, original=previous[case], fresh=fresh))
        comparisons.append(dict(comparison, independent_frequency_error_cm=frequency_error,
            independent_nuclear_error_eV=nuclear_error, independent_G_error_eV=abs(G-status['G_eV'])))
    if {r['case'] for r in states} != {'formaldehyde', 'methyl_radical'}:
        raise ValueError('Both predeclared spin controls are required')
    out.mkdir(parents=True, exist_ok=False)
    atomic_write_json(out / 'states.json', dict(states=states))
    result = dict(status='fresh_compatible_GPU_deployment_independently_checked', checked=time.time(),
        controls=comparisons, tolerances=report['tolerances'], source_commit=installation['source_commit'],
        external_distributions=installation['third_party_distributions'],
        recorded_environment_dependency_conflicts=14, compatible_CUDA_runtime_changes=14,
        libxtb_sha256=installation['libxtb_sha256'], original_environment_modified=False,
        method_packages_changed=False, full_portability_claimed=False,
        source_sha256={str(p):digest(p) for p in (Path(__file__).resolve(), campaign/'summary.json',
            campaign/'plan.json', restoration/'status.json', evidence/'reference_control_snapshot.json',
            evidence/'certificate_closure.json', out/'states.json')},
        original_and_fresh_backing_sha256={**reference['source_sha256'], **report['backing_sha256']},
        limits=report['limits'])
    atomic_write_json(out / 'summary.json', result)
    print(json.dumps(dict(status=result['status'], controls=len(states))))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('restoration', 'campaign', 'evidence', 'out'):
        parser.add_argument('--'+name, type=Path, required=True)
    args = parser.parse_args()
    publish(args.restoration.resolve(), args.campaign.resolve(), args.evidence.resolve(), args.out.resolve())


if __name__ == '__main__':
    main()
