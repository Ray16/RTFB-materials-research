#!/usr/bin/env python
"""Retain a failed native Hessian and compare its worst asymmetry to gradients.

Read-only with respect to production. No Hessian is symmetrized or admitted here.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import socket
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'src'))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--manifest', type=Path, required=True)
    ap.add_argument('--task', type=int, required=True)
    a = ap.parse_args()
    from redox.validation.spins import atomic_json
    from redox.qm.native_thermal import NATIVE_THERMAL, implementation_controls
    from redox.qm.dft import _build_smd_mf, _ecp_for, _orbital_diag, _free_gpu_pool
    from ase.io import read
    from pyscf import gto
    import numpy as np
    import gpu4pyscf.df.hessian.rhf as dr
    import gpu4pyscf.df.hessian.uhf as du
    manifest = json.loads(a.manifest.read_text())
    for path, digest in manifest['producer_sha256'].items():
        assert hashlib.sha256(Path(path).read_bytes()).hexdigest() == digest
    task = manifest['tasks'][a.task]
    assert task['kind'] == 'hessian_diagnostic'
    assert os.getenv('CUDA_VISIBLE_DEVICES') and sorted(os.sched_getaffinity(0)) == task['compute_cpus']
    out = a.manifest.parent / task['id']
    out.mkdir()
    source = Path(task['source'])
    raw = json.loads((source / 'result.json').read_text())
    atoms = read(source / 'opt.xyz')
    symbols = atoms.get_chemical_symbols()
    coords = atoms.positions / .52917721067
    sha = hashlib.sha256((source / 'opt.xyz').read_bytes()).hexdigest()
    (out / 'source.xyz').write_bytes((source / 'opt.xyz').read_bytes())
    report = dict(status='running', task=task, host=socket.gethostname(), pid=os.getpid(),
                  started=time.time(), protocol=NATIVE_THERMAL, source_sha256=sha,
                  controls=implementation_controls(), production_modified=False)
    def save():
        atomic_json(out / 'status.json', report)
    save()
    try:
        for module in (dr, du):
            if hasattr(module, 'get_avail_mem'):
                old = module.get_avail_mem
                module.get_avail_mem = lambda *args, _old=old, **kwargs: min(_old(*args, **kwargs), 6 * 1024**3)
        def solve(x):
            mol = gto.M(atom=list(zip(symbols, x.tolist())), unit='Bohr', charge=raw['charge'],
                        spin=raw['mult'] - 1, basis='def2-svpd' if raw['charge'] < 0 else 'def2-svp',
                        ecp=_ecp_for(symbols, 'def2-svpd' if raw['charge'] < 0 else 'def2-svp'),
                        verbose=0, max_memory=6000)
            mf = _build_smd_mf(mol, 'r2scan', 'acetonitrile', 'gpu', disp='d4')
            mf.grids.level = 3
            mf.conv_tol = 1e-12
            mf.conv_tol_grad = 1e-8
            mf.max_cycle = 200
            mf.kernel()
            assert mf.converged
            return mf
        def cpu(x):
            return np.asarray(x.get() if hasattr(x, 'get') else x)
        mf = solve(coords)
        report['spin_diagnostics'] = _orbital_diag(mf, raw['mult'] - 1)
        grad = mf.nuc_grad_method()
        grad.grid_response = True
        report['gradient_max'] = float(abs(cpu(grad.kernel())).max())
        mf.cphf_grids = mf.grids.copy()
        ho = mf.Hessian()
        ho.grid_response = True
        ho.max_cycle = 100
        report['status'] = 'computing_full_response_hessian'
        save()
        H = cpu(ho.kernel())
        mat = H.transpose(0, 2, 1, 3).reshape(len(coords)*3, len(coords)*3)
        np.save(out / 'hessian_unsymmetrized.npy', mat)
        asym = abs(mat - mat.T)
        i, j = np.unravel_index(np.nanargmax(asym), asym.shape)
        report.update(maximum_asymmetry=float(asym[i, j]), worst_indices=[int(i), int(j)],
                      all_finite=bool(np.isfinite(mat).all()), analytic_pair=[float(mat[i,j]), float(mat[j,i])],
                      status='checking_conservative_gradients', checks=[])
        save()
        del mf, ho, grad
        _free_gpu_pool()
        for column in sorted({int(i), int(j)}):
            for step in (.0005, .001):
                gs = []
                for sign in (-1, 1):
                    x = coords.copy()
                    x.flat[column] += sign * step
                    fm = solve(x)
                    gm = fm.nuc_grad_method()
                    gm.grid_response = True
                    gs.append(cpu(gm.kernel()).reshape(-1))
                    del fm, gm
                    _free_gpu_pool()
                numerical = (gs[1] - gs[0]) / (2*step)
                # Native storage uses [perturbed coordinate, gradient coordinate].
                error = float(np.linalg.norm(numerical - mat[column,:]) / np.linalg.norm(numerical))
                report['checks'].append(dict(column=column, step_Bohr=step, relative_error=error,
                                             numerical_pair=float(numerical[j if column == i else i])))
                save()
        report.update(status='diagnostic_complete', source_unchanged=sha == hashlib.sha256((source/'opt.xyz').read_bytes()).hexdigest(),
                      finished=time.time(), production_admission=False)
        save()
    except BaseException as exc:
        report.update(status='failed_preserved', error=repr(exc), finished=time.time())
        save()
        raise


if __name__ == '__main__':
    main()
