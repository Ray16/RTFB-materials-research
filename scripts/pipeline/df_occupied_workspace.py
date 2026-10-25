"""Bound occupied-orbital blocks in the checked DF derivative execution.

Caps 4 and 32 exercise the same exact installed loop. Full replay at both caps
is required; the smaller cap forces multiple blocks in the singlet control.
"""
from contextlib import contextmanager
import hashlib
import inspect
from pathlib import Path
import textwrap
import time


def cap_source(body):
    changes = (
        ('    nocc_per_batch = min(nocc, (nocc_per_batch + 7) & 0xffffff8)',
         '    nocc_per_batch = min(nocc, _occupied_cap, (nocc_per_batch + 7) & 0xffffff8)'),
        ('            nocc_in_batch = nocc1 - nocc0',
         "            nocc_in_batch = nocc1 - nocc0\n"
         "            _occupied_trace(dict(batch=int(batch_id), begin=int(nocc0), end=int(nocc1),\n"
         "                size=int(nocc_in_batch), total=int(nocc)))"),
    )
    for old, new in changes:
        if body.count(old) != 1:
            raise ValueError('Installed occupied-orbital batching template changed')
        body = body.replace(old, new, 1)
    return body


@contextmanager
def compact_smd_execution(directory, budget_GiB, occupied_cap=32):
    if type(occupied_cap) is not int or occupied_cap not in (4, 32):
        raise ValueError('Only the two replayed occupied caps are supported')
    import cupy as cp
    import numpy as np
    from gpu4pyscf.df.hessian import rhf, rks, uks
    from redox.core.common import ROOT
    from redox.core.io import atomic_write_json
    from redox.workflows import polymer_campaign as api

    directory = Path(directory)
    parent_path = ROOT / 'scripts/pipeline/df_derivative_workspace.py'
    parent = api.load_module('occupied_block_DF_underlay', parent_path)
    stats = dict(status='occupied_block_execution_active_unqualified', started=time.time(),
        occupied_cap=occupied_cap, batches=[], checkpoints=[], precision='float64', equations_changed=False,
        source_sha256={str(p): api.digest(p) for p in (Path(__file__).resolve(), parent_path, Path(rhf.__file__))})
    with parent.compact_smd_execution(directory, budget_GiB) as proof:
        old = rhf._get_veff
        clone = getattr(old, '__kwdefaults__', {}).get('_function')
        if clone is None:
            raise ValueError('The checked derivative profile must expose its exact function')
        body = cap_source(inspect.getsource(clone))
        path = directory / 'occupied_capped_get_veff.py'
        path.write_text(body)
        def trace(row):
            stats['batches'].append(dict(checked=time.time(), **row))
            atomic_write_json(directory / 'occupied_progress.json', stats)
        namespace = dict(clone.__globals__, _occupied_cap=occupied_cap, _occupied_trace=trace)
        local = {}
        exec(compile(textwrap.dedent(body), str(path), 'exec'), namespace, local)
        modified = local['_get_veff']
        def wrapper(*args, **kwargs):
            return old(*args, _function=modified, **kwargs)
        rhf._get_veff = wrapper
        originals = []
        def digest_array(value):
            array = cp.asnumpy(value) if isinstance(value, cp.ndarray) else np.asarray(value)
            return dict(shape=list(array.shape), dtype=str(array.dtype), sha256=hashlib.sha256(array.tobytes()).hexdigest())
        for spin, module in (('rks', rks), ('uks', uks)):
            for name in ('partial_hess_elec', 'make_h1'):
                original = getattr(module.Hessian, name); tag = spin + '_' + name
                def checkpoint(hess, *args, _original=original, _tag=tag, **kwargs):
                    mf = hess.base
                    identity = dict(stage=_tag, coordinates_Bohr=mf.mol.atom_coords().tolist(),
                        basis=repr(mf.mol.basis), charge=mf.mol.charge, mult=mf.mol.spin + 1,
                        xc=mf.xc, SCF_E_Ha=float(mf.e_tot), SCF_converged=bool(mf.converged),
                        mo_coeff=digest_array(mf.mo_coeff), mo_energy=digest_array(mf.mo_energy),
                        mo_occ=digest_array(mf.mo_occ), qualified_native_G=False)
                    value = _original(hess, *args, **kwargs)
                    identity['arrays'] = {str(p): api.digest(p) for p in directory.glob(_tag + '_*.npy')}
                    target = directory / (_tag + '_state_identity.json')
                    atomic_write_json(target, identity)
                    stats['checkpoints'].append(dict(stage=_tag, path=str(target), sha256=api.digest(target)))
                    return value
                originals.append((module.Hessian, name, original))
                setattr(module.Hessian, name, checkpoint)
        proof['source_sha256'].update(stats['source_sha256'])
        proof['DF_derivative_occupied_cap'] = occupied_cap
        try:
            yield proof
            stats['status'] = 'occupied_block_execution_finished_requires_replay_and_physical_checks'
        except BaseException as error:
            stats.update(status='occupied_block_execution_failed_preserved', error=repr(error))
            raise
        finally:
            for cls, name, original in originals:
                setattr(cls, name, original)
            rhf._get_veff = old
            stats['finished'] = time.time()
            atomic_write_json(directory / 'occupied_workspace.json', stats)
