#!/usr/bin/env python
"""Compare complete engine replays to the preserved einsum reference."""
import argparse
import json
from pathlib import Path
import time

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api


def checked_controls(root):
    path = Path(root) / 'comparison.json'; row = json.loads(path.read_text())
    if (row['status'] != 'two_spin_full_cutensor_and_fallback_controls_match'
            or row['installed_environment_modified'] or row['large_native_G_qualified']
            or len(row['controls']) != 2):
        raise ValueError('Both complete cross-engine controls must pass within scope')
    api.verify_sources(row['source_sha256'])
    parent = api.load_module('cutensor_retained_occupied_cap_proofs', ROOT / 'scripts/validation/qualify_df_occupied_workspace.py')
    proof = parent.checked_controls(Path(row['fallback_root']))
    proof.update(row['source_sha256']); proof[str(path)] = api.digest(path)
    return proof


def review(root, fallback, profile):
    import numpy as np
    definition = json.loads(profile.read_text())
    parent = api.load_module('cutensor_fallback_full_cap_checks', ROOT / 'scripts/validation/qualify_df_occupied_workspace.py')
    proof = parent.checked_controls(fallback)
    proof.update({str(p): api.digest(p) for p in (Path(__file__).resolve(), profile,
        ROOT / 'scripts/fleet/tensor_runtime_worker.py')})
    rows = []
    for mult in (1, 2):
        folder = root / 'caps/32/controls' / str(mult)
        old_folder = fallback / 'caps/32/controls' / str(mult)
        current_path, old_path = folder / 'status.json', old_folder / 'status.json'
        current, old = (json.loads(p.read_text()) for p in (current_path, old_path))
        if current['status'] != 'full_original_compact_equivalence_passed' or not current.get('occupied_block_checks_passed'):
            raise ValueError('New full replay has not passed')
        api.verify_sources(current['source_sha256']); api.verify_sources(current['evidence_sha256'])
        runtime_path = root / 'runtime' / str(mult) / 'status.json'
        runtime = json.loads(runtime_path.read_text())
        if (runtime['status'] != 'task_returned_under_checked_cutensor_runtime' or runtime['engine'] != 'cutensor'
                or runtime['versions'] != definition['versions'] or runtime['library_sha256'] != definition['library_sha256']
                or runtime['profile_sha256'] != api.digest(profile)
                or runtime['source_sha256'] != api.digest(ROOT / 'scripts/fleet/tensor_runtime_worker.py')
                or runtime['pid'] != current['pid']):
            raise ValueError('Actual contraction runtime is not bound to this full replay')
        guard = root / 'guards' / str(mult)
        owner, state = (json.loads((guard / n).read_text()) for n in ('owner.json', 'state.json'))
        if (owner['pid'] != current['pid'] or not owner.get('guard_stopped')
                or not state['own_residency_observed'] or state['foreign_resident_PIDs']
                or (guard / 'interruption.json').exists()):
            raise ValueError('Full replay residency was not exclusive')
        differences = {}
        for arm in ('original', 'compact'):
            values, reference = current['arms'][arm], old['arms']['original']
            differences[arm] = dict(
                energy_eV=max(abs(values[k] - reference[k]) for k in ('E_smd_eV', 'E_gas_eV')),
                geometry_level_energy_Ha=abs(values['geometry_level_E_Ha'] - reference['geometry_level_E_Ha']),
                coordinates_A=float(np.max(abs(np.array(values['coordinates_A']) - np.array(reference['coordinates_A'])))),
                gradient=float(np.max(abs(np.load(folder / arm / 'gradient.npy') - np.load(old_folder / 'original/gradient.npy')))),
                full_Hessian=float(np.max(abs(np.load(folder / arm / 'hessian.npy') - np.load(old_folder / 'original/hessian.npy')))))
            if any(not np.isfinite(v) or v > 1e-6 for v in differences[arm].values()):
                raise ValueError('Full cross-engine replay exceeds the original 1e-6 bounds')
        for p in (current_path, old_path, runtime_path, guard / 'owner.json', guard / 'state.json'):
            proof[str(p)] = api.digest(p)
        proof.update(current['source_sha256']); proof.update(current['evidence_sha256'])
        rows.append(dict(mult=mult, differences=differences, actual_engine='cutensor',
            peak_GPU_MiB=current['observed_peak_GPU_MiB'], device=runtime['device_name']))
    result = dict(status='two_spin_full_cutensor_and_fallback_controls_match', checked=time.time(),
        controls=rows, fallback_root=str(fallback), runtime_profile=str(profile), source_sha256=proof,
        installed_environment_modified=False, method_packages_changed=False,
        large_native_G_qualified=False, production_default_changed=False)
    atomic_write_json(root / 'comparison.json', result)
    print(json.dumps(dict(status=result['status'], controls=rows)))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('root', 'fallback', 'profile'):
        parser.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args()
    from redox.workflows.resources import allocation_guard, cpu_allocation, host_capacity
    import os
    with allocation_guard():
        host_capacity(1); os.sched_setaffinity(0, cpu_allocation(1))
    review(args.root.resolve(), args.fallback.resolve(), args.profile.resolve())


if __name__ == '__main__':
    main()
