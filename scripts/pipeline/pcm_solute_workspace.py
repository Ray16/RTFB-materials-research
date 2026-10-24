"""Suspend unused solvent matrices during pure solute Hessian subroutines.

The full solvent matrices return byte-exactly before solvent response. Access
during suspension fails closed. No integral, derivative, grid or precision is
changed. Full singlet/doublet Hessian replay is required before production.
"""
from contextlib import contextmanager
import gc
from pathlib import Path
import time

import numpy as np


class SuspendedIntermediates(dict):
    def __init__(self, available, suspended):
        super().__init__(available)
        self.suspended = frozenset(suspended)

    def __getitem__(self, key):
        if key in self.suspended:
            raise RuntimeError('Pure solute subroutine accessed suspended solvent matrix: ' + key)
        return super().__getitem__(key)

    def get(self, key, default=None):
        if key in self.suspended:
            raise RuntimeError('Pure solute subroutine accessed suspended solvent matrix: ' + key)
        return super().get(key, default)


@contextmanager
def suspend_matrices(solvent, xp, metrics, release):
    original = solvent._intermediates
    keys = ('S', 'D', 'K_LU')
    if original is None or any(key not in original for key in keys):
        raise ValueError('A built full SMD surface is required')
    if any(original[key].dtype != xp.float64 or original[key].ndim != 2 for key in keys):
        raise ValueError('Expected FP64 solvent matrices')
    host = {}
    metrics['suspensions'] += 1
    try:
        for key in keys:
            value = original[key]
            host[key] = xp.asnumpy(value)
            metrics['maximum_suspended_bytes'] = max(metrics['maximum_suspended_bytes'], sum(v.nbytes for v in host.values()))
            del original[key], value
        unavailable = SuspendedIntermediates(original, keys)
        solvent._intermediates = unavailable
        release()
        yield
        if solvent._intermediates is not unavailable:
            raise RuntimeError('Solvent intermediates changed during pure solute Hessian work')
    finally:
        release()
        # Exact value checks qualify storage changes, including restored LU
        # factors; no factorization or density solve is substituted.
        for key in tuple(host):
            value = xp.asarray(host[key])
            restored = xp.asnumpy(value)
            if (restored.shape != host[key].shape or restored.dtype != host[key].dtype
                    or not np.array_equal(np.ascontiguousarray(restored).view(np.uint8),
                                          np.ascontiguousarray(host[key]).view(np.uint8))):
                raise ValueError('Solvent matrix round trip changed values: ' + key)
            original[key] = value
            del host[key], value, restored
            metrics['exact_matrix_round_trips'] += 1
        solvent._intermediates = original
        metrics['restorations'] += 1


@contextmanager
def compact_smd_execution(directory, budget_GiB):
    import cupy as cp
    from gpu4pyscf.df.hessian import rks, uks
    from gpu4pyscf.solvent import smd
    from redox.core.common import ROOT
    from redox.core.io import atomic_write_json
    from redox.workflows import polymer_campaign as api

    directory = Path(directory)
    underlay = api.load_module('solute_workspace_streamed_weighted_parent',
        ROOT / 'scripts/pipeline/streamed_pcm_weighted.py')
    stats = dict(status='solute_workspace_active_unqualified', started=time.time(),
        suspensions=0, restorations=0, exact_matrix_round_trips=0, maximum_suspended_bytes=0,
        calls={label + '.' + name: 0 for label in ('rks', 'uks') for name in ('partial_hess_elec', 'make_h1')},
        equations_changed=False, precision_changed=False,
        source_sha256={str(p): api.digest(p) for p in (Path(__file__).resolve(), Path(rks.__file__), Path(uks.__file__))})

    def release():
        gc.collect()
        cp.get_default_memory_pool().free_all_blocks()

    with underlay.compact_smd_execution(directory, budget_GiB) as proof:
        originals = []
        for label, module in (('rks', rks), ('uks', uks)):
            for name in ('partial_hess_elec', 'make_h1'):
                original = getattr(module.Hessian, name)
                tag = label + '.' + name
                def wrapped(hess, *args, _old=original, _tag=tag, **kwargs):
                    solvent = getattr(hess.base, 'with_solvent', None)
                    if not isinstance(solvent, smd.SMD) or solvent.method.upper() != 'SMD':
                        return _old(hess, *args, **kwargs)
                    stats['calls'][_tag] += 1
                    with suspend_matrices(solvent, cp, stats, release):
                        return _old(hess, *args, **kwargs)
                originals.append((module.Hessian, name, original))
                setattr(module.Hessian, name, wrapped)
        proof['source_sha256'].update(stats['source_sha256'])
        proof.update(solute_workspace_solvent_suspension=True, solvent_response_omitted=False)
        try:
            yield proof
            stats['status'] = 'solute_workspace_finished_requires_full_controls_and_physical_review'
        except BaseException as error:
            stats.update(status='solute_workspace_failed_preserved', error=repr(error))
            raise
        finally:
            for cls, name, original in originals:
                setattr(cls, name, original)
            stats['finished'] = time.time()
            atomic_write_json(directory / 'solute_workspace.json', stats)
