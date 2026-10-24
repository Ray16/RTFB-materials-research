"""Bound CPU PCM integral tiles through the existing solver's memory setting.

PySCF retains a 400-point minimum tile. The molecule/DF memory budget and
solvent tile budget are deliberately separate. All integral and contraction
functions remain unchanged; wrappers only set the budget and record allocation.
"""
from contextlib import contextmanager
import inspect
import math
from pathlib import Path
import time

from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api


@contextmanager
def bounded_pcm_memory(directory, maximum_PCM_memory_MiB=8192):
    if not math.isfinite(maximum_PCM_memory_MiB) or maximum_PCM_memory_MiB <= 0:
        raise ValueError('Positive independent PCM memory budget required')
    from pyscf import df
    from pyscf.solvent import pcm
    from redox.qm import dft
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=False)
    original_build = dft._build_smd_mf
    original_v, original_vmat = pcm.PCM._get_v, pcm.PCM._get_vmat
    original_aux = df.incore.aux_e2
    active = []
    report = dict(status='CPU_PCM_budget_active', started=time.time(),
        PCM_memory_MiB=maximum_PCM_memory_MiB, molecule_DF_budget_changed=False,
        integral_or_contraction_equations_changed=False, minimum_tile_points=400,
        builds=[], tiles={'potential':0,'matrix':0}, maximum_tile_points=0,
        maximum_integral_tensor_bytes=0, source_sha256={str(p):api.digest(p) for p in
            (Path(__file__).resolve(), Path(pcm.__file__), Path(df.incore.__file__), Path(dft.__file__))})

    def build(*a, **kw):
        mf = original_build(*a, **kw)
        previous = mf.with_solvent.max_memory
        mf.with_solvent.max_memory = min(previous, maximum_PCM_memory_MiB)
        report['builds'].append(dict(atoms=mf.mol.natm, nao=mf.mol.nao_nr(),
            original_PCM_memory_MiB=previous, PCM_memory_MiB=mf.with_solvent.max_memory,
            molecule_memory_MiB=mf.mol.max_memory))
        atomic_write_json(directory/'status.json', report)
        return mf

    def potential(self, *a, **kw):
        active.append('potential')
        try:
            return original_v(self, *a, **kw)
        finally:
            active.pop()

    def matrix(self, *a, **kw):
        active.append('matrix')
        try:
            return original_vmat(self, *a, **kw)
        finally:
            active.pop()

    def integrals(*a, **kw):
        value = original_aux(*a, **kw)
        if active:
            report['tiles'][active[-1]] += 1
            report['maximum_tile_points'] = max(report['maximum_tile_points'], value.shape[-1])
            report['maximum_integral_tensor_bytes'] = max(report['maximum_integral_tensor_bytes'], value.nbytes)
            atomic_write_json(directory/'status.json', report)
        return value

    for name, function in (('get_v',original_v),('get_vmat',original_vmat)):
        (directory/(name+'_unmodified_source.py')).write_text(inspect.getsource(function))
    dft._build_smd_mf = build
    pcm.PCM._get_v, pcm.PCM._get_vmat = potential, matrix
    df.incore.aux_e2 = integrals
    try:
        yield report
        report['status'] = 'CPU_PCM_budget_finished_requires_numerical_review'
    except BaseException as error:
        report.update(status='CPU_PCM_budget_failed_preserved', error=repr(error))
        raise
    finally:
        dft._build_smd_mf = original_build
        pcm.PCM._get_v, pcm.PCM._get_vmat = original_v, original_vmat
        df.incore.aux_e2 = original_aux
        report['finished'] = time.time()
        atomic_write_json(directory/'status.json', report)
