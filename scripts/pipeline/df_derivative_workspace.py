"""Use the installed streamed AO-derivative path with bounded DF batches.

Full spin/Hessian replay and a large-state resource check are mandatory. This
selects an existing contraction path, retaining every integral and FP64 value.
"""
from contextlib import contextmanager
import inspect
from pathlib import Path
import textwrap
import time


def transformed_source(body, spin):
    if spin not in ('rhf', 'uhf'):
        raise ValueError('Only the installed restricted/unrestricted DF kernels are supported')
    expression = 'natm*3*nao*nao' if spin == 'rhf' else 'natm*3*2*nao**2'
    changes = (
        (f'    mem_sufficient = {expression} < word_avail * .6',
         '    mem_sufficient = False # Retain the installed MO-contracted derivative path'),
        ('    batch_size = int(word_avail*.8 / _unit)',
         '    batch_size = int(word_avail*.25 / _unit)'),
        ("    assert blksize > 0, 'Insufficient GPU memory'",
         "    _workspace_trace(dict(mem_avail=int(mem_avail), word_avail=int(word_avail),\n"
         "        natm=int(natm), nao=int(nao), nocc=int(nocc), naux=int(naux),\n"
         "        nao_pair=int(nao_pair), auxiliary_batch=int(batch_size), blksize=int(blksize),\n"
         "        full_atom_AO_cache=bool(mem_sufficient)))\n"
         "    assert blksize > 0, 'Insufficient GPU memory'"),
    )
    for old, new in changes:
        if body.count(old) != 1:
            raise ValueError('Installed DF derivative allocation template changed')
        body = body.replace(old, new, 1)
    return body


@contextmanager
def compact_smd_execution(directory, budget_GiB):
    import cupy as cp
    import numpy as np
    from gpu4pyscf.df.hessian import rhf, uhf, rks, uks
    from redox.core.common import ROOT
    from redox.core.io import atomic_write_json
    from redox.workflows import polymer_campaign as api

    directory = Path(directory)
    parent_path = ROOT / 'scripts/pipeline/pcm_solute_workspace.py'
    parent = api.load_module('DF_derivative_solute_workspace_parent', parent_path)
    stats = dict(status='DF_derivative_workspace_active_unqualified', started=time.time(),
        calls={'rhf': 0, 'uhf': 0}, allocations=[], checkpoints=[], equations_changed=False,
        precision_changed=False, full_atom_AO_cache=False, auxiliary_budget_fraction=.25,
        source_sha256={str(p): api.digest(p) for p in
            (Path(__file__).resolve(), parent_path, Path(rhf.__file__), Path(uhf.__file__))})
    with parent.compact_smd_execution(directory, budget_GiB) as proof:
        functions, methods = [], []
        def trace(spin, row):
            stats['allocations'].append(dict(spin=spin, checked=time.time(), **row))
            atomic_write_json(directory / 'DF_derivative_progress.json', stats)
        for spin, module in (('rhf', rhf), ('uhf', uhf)):
            original = module._get_veff
            body = transformed_source(inspect.getsource(original), spin)
            path = directory / (spin + '_bounded_get_veff.py')
            header = '\n'.join(Path(module.__file__).read_text().splitlines()[:13]) + '\n\n'
            path.write_text(header + body)
            scope = dict(original.__globals__, _workspace_trace=lambda row, _spin=spin: trace(_spin, row))
            local = {}
            exec(compile(textwrap.dedent(body), str(path), 'exec'), scope, local)
            replacement = local['_get_veff']
            def wrapped(*args, _function=replacement, _spin=spin, **kwargs):
                stats['calls'][_spin] += 1
                return _function(*args, **kwargs)
            functions.append((module, original))
            module._get_veff = wrapped
        for spin, module in (('rks', rks), ('uks', uks)):
            for name in ('partial_hess_elec', 'make_h1'):
                original = getattr(module.Hessian, name)
                tag = spin + '_' + name
                def checkpoint(hess, *args, _function=original, _tag=tag, **kwargs):
                    atomic_write_json(directory / 'DF_stage.json', dict(stage=_tag, status='running',
                        updated=time.time(), qualified_native_G=False))
                    result = _function(hess, *args, **kwargs)
                    parts = result if isinstance(result, (tuple, list)) else [result]
                    for index, part in enumerate(parts):
                        target = directory / (_tag + '_' + str(index) + '.npy')
                        value = cp.asnumpy(part) if isinstance(part, cp.ndarray) else np.asarray(part)
                        np.save(target, value)
                        item = dict(stage=_tag, part=index, path=str(target), sha256=api.digest(target),
                            shape=list(value.shape), coordinates_Bohr=hess.mol.atom_coords().tolist(),
                            charge=hess.mol.charge, mult=hess.mol.spin + 1, qualified_native_G=False)
                        stats['checkpoints'].append(item)
                        del value
                    atomic_write_json(directory / 'DF_stage.json', dict(stage=_tag,
                        status='checkpointed', parts=len(parts), updated=time.time(), qualified_native_G=False))
                    return result
                methods.append((module.Hessian, name, original))
                setattr(module.Hessian, name, checkpoint)
        proof['source_sha256'].update(stats['source_sha256'])
        proof.update(DF_derivative_full_atom_AO_cache=False, DF_derivative_auxiliary_fraction=.25)
        try:
            yield proof
            stats['status'] = 'DF_derivative_workspace_finished_requires_full_replay_and_physical_checks'
        except BaseException as error:
            stats.update(status='DF_derivative_workspace_failed_preserved', error=repr(error))
            raise
        finally:
            for cls, name, original in methods:
                setattr(cls, name, original)
            for module, original in functions:
                module._get_veff = original
            stats['finished'] = time.time()
            atomic_write_json(directory / 'DF_derivative_workspace.json', stats)
