"""CPU-only execution settings and durable SCF records; quantum models unchanged."""
from contextlib import contextmanager
from importlib import metadata
import json
from pathlib import Path
import sys
import time

from redox.core.io import atomic_write_json


@contextmanager
def cpu_execution(directory, memory_MiB):
    import pyscf
    from pyscf import gto
    from redox.core import protocol
    from redox.qm import dft
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=False)
    original_M, original_kernel = gto.M, dft._kernel_robust
    original_free, original_versions = dft._free_gpu_pool, protocol.software_versions
    report = dict(started=time.time(), status='CPU_execution_active', maximum_library_memory_MiB=memory_MiB,
        scientific_model_changed=False, GPU_imports_requested=False, SCF_calls=[])

    def molecule(*args, **kwargs):
        kwargs['max_memory'] = memory_MiB
        return original_M(*args, **kwargs)

    def versions():
        # Distribution metadata does not initialize the installed CUDA backend.
        return dict(pyscf=metadata.version('pyscf'),
            gpu4pyscf=metadata.version('gpu4pyscf-cuda12x'))

    def kernel(mf, dm0=None):
        index = len(report['SCF_calls'])
        folder = directory / f'SCF_{index:03d}'
        folder.mkdir()
        mf.chkfile = str(folder / 'SCF.h5')
        progress = dict(qualified_result=False, SCF_converged=False, xc=mf.xc,
            basis=repr(mf.mol.basis), charge=mf.mol.charge, mult=mf.mol.spin+1,
            solvent=getattr(getattr(mf, 'with_solvent', None), 'solvent', None), started=time.time())
        before = getattr(mf, 'callback', None)
        def callback(environment):
            if callable(before):
                before(environment)
            progress.update(cycle=int(environment['cycle'])+1, E_Ha=float(environment['e_tot']), updated=time.time())
            atomic_write_json(folder / 'progress.json', progress)
        mf.callback = callback
        report['SCF_calls'].append(str(folder))
        atomic_write_json(directory / 'status.json', report)
        try:
            value = original_kernel(mf, dm0=dm0)
            progress.update(SCF_converged=bool(value[1]), E_Ha=float(value[0]), finished=time.time())
            return value
        finally:
            mf.callback = before
            atomic_write_json(folder / 'progress.json', progress)

    gto.M, dft._kernel_robust = molecule, kernel
    dft._free_gpu_pool = lambda: None
    protocol.software_versions = versions
    try:
        yield report
        if any(name == 'cupy' or name.startswith('cupy.') or name == 'gpu4pyscf'
                or name.startswith('gpu4pyscf.') for name in sys.modules):
            raise RuntimeError('CPU execution imported a CUDA backend')
        report['status'] = 'CPU_execution_finished_requires_physical_review'
    except BaseException as error:
        report.update(status='CPU_execution_failed_preserved', error=repr(error))
        raise
    finally:
        gto.M, dft._kernel_robust = original_M, original_kernel
        dft._free_gpu_pool, protocol.software_versions = original_free, original_versions
        report['finished'] = time.time()
        atomic_write_json(directory / 'status.json', report)
