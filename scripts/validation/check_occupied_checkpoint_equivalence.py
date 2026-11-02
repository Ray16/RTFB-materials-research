#!/usr/bin/env python3
"""Check compact storage against actual pre-existing full-orbital determinants.

This CPU-only control does not rerun quantum calculations or qualify new G.
"""
import argparse
import json
import os
from pathlib import Path

import numpy as np
from pyscf.lib import chkfile
from pyscf.scf.hf import make_rdm1

from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json
from redox.qm.occupied_checkpoints import save, read, density, review
from redox.validation.scf_frame import checkpoint_spin


def check(sources, output):
    if output.exists():
        raise FileExistsError('Preserve existing equivalence controls')
    output.mkdir(parents=True); records = []
    for i, path in enumerate(sources):
        mol = chkfile.load_mol(str(path)); full = chkfile.load(str(path), 'scf')
        target = output / f'{i:03d}.h5'; original_spin = checkpoint_spin(path)
        save(target, mol, full['mo_coeff'], full['mo_occ'], full['mo_energy'], full['e_tot'],
             converged=False, cycle=-1)
        _, stored = read(target); proof = review(target)
        identical = np.array_equal(stored['mo_coeff_occupied'], full['mo_coeff'][:, full['mo_occ'] > 0])
        density_error = float(np.max(abs(density(target) - make_rdm1(full['mo_coeff'], full['mo_occ']))))
        if (not identical or not np.array_equal(stored['mo_occ'], full['mo_occ'])
                or not np.array_equal(stored['mo_energy'], full['mo_energy'])
                or stored['e_tot'] != full['e_tot'] or density_error > 1e-13
                or abs(original_spin['s_squared'] - proof['restricted_s_squared']) > 1e-7):
            raise ValueError('Compact storage changed an actual source determinant/density/spin')
        records.append(dict(source=str(path), source_sha256=sha256_file(path),
            occupied_checkpoint=str(target), occupied_checkpoint_sha256=sha256_file(target),
            actual_molecular_frame_preserved=True, occupied_coefficients_bitwise_identical=identical,
            all_occupations_and_orbital_energies_bitwise_identical=True,
            energy_bitwise_identical=True, maximum_density_error=density_error, **proof))
    import redox.qm.occupied_checkpoints as implementation
    result = dict(status='actual_saved_determinants_preserved_by_occupied_storage', records=records,
        source_sha256={str(p): sha256_file(p) for p in [Path(__file__), Path(implementation.__file__)]},
        new_SCF_calculations=False, source_convergence_requalified=False,
        native_G_qualified=False, full_response_restart_supported=False)
    atomic_write_json(output / 'report.json', result)
    return result


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source', type=Path, action='append', required=True)
    p.add_argument('--out', type=Path, required=True); a = p.parse_args()
    from redox.workflows.resources import allocation_guard, host_capacity, cpu_allocation
    with allocation_guard():
        host_capacity(1, minimum_ram_gib=4); os.sched_setaffinity(0, cpu_allocation(1))
    result = check([p.resolve() for p in a.source], a.out.resolve())
    print(json.dumps(dict(status=result['status'], actual_determinants=len(result['records']),
        maximum_density_error=max(r['maximum_density_error'] for r in result['records']))))
