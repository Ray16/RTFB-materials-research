#!/usr/bin/env python
"""Fail closed on a declared contraction engine before running an owned task."""
import argparse
import hashlib
import importlib.metadata as metadata
import json
import os
from pathlib import Path
import runpy
import sys
import time

from redox.core.io import atomic_write_json
from redox.workflows.resources import verify_gate_owner


def digest(path):
    result = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024**2), b''):
            result.update(block)
    return result.hexdigest()


def validate(profile):
    import cupy as cp
    import numpy as np
    import gpu4pyscf.lib.cutensor as engine
    if engine.cutensor is None or engine.contract_engine is not None:
        raise RuntimeError('The declared cuTENSOR engine is unavailable; refuse silent einsum fallback')
    versions = {name: metadata.version(name) for name in profile['versions']}
    if versions != profile['versions']:
        raise ValueError('Contraction runtime distributions differ from the qualified profile')
    library = Path(metadata.distribution('cutensor-cu12').locate_file('cutensor/lib'))
    for name, expected in profile['library_sha256'].items():
        if digest(library / name) != expected:
            raise ValueError('Contraction library differs: ' + name)
    rng = np.random.default_rng(20261025)
    cases = [('ij,jk->ik', rng.normal(size=(5, 7)), rng.normal(size=(7, 4))),
        ('xrs,yts->rtxy', rng.normal(size=(3, 5, 7)), rng.normal(size=(3, 4, 7))),
        ('uv,uij->vij', rng.normal(size=(7, 5)), rng.normal(size=(7, 3, 4)))]
    errors = []
    for pattern, a, b in cases:
        reference = np.einsum(pattern, a, b)
        out = cp.ones(reference.shape, dtype=cp.float64)
        got = engine.contract(pattern, cp.asarray(a), cp.asarray(b), alpha=1.5, beta=.25, out=out)
        error = float(np.max(abs(cp.asnumpy(got) - (1.5 * reference + .25))))
        if got is not out or got.dtype != cp.float64 or not np.isfinite(error) or error > 1e-11:
            raise ValueError('FP64 contraction/out/alpha/beta smoke check failed')
        errors.append(dict(pattern=pattern, maximum_absolute_error=error))
    device = cp.cuda.runtime.getDeviceProperties(0)
    return dict(engine='cutensor', versions=versions, FP64_smoke_checks=errors,
        library_sha256=profile['library_sha256'], actual_library_directory=str(library),
        device_name=device['name'].decode() if isinstance(device['name'], bytes) else str(device['name']),
        compute_capability=[device['major'], device['minor']],
        total_GPU_bytes=int(device['totalGlobalMem']), source_sha256=digest(__file__),
        installed_packages_changed=False, validated_full_molecular_workflow=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--profile', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('program', type=Path)
    parser.add_argument('arguments', nargs=argparse.REMAINDER)
    args = parser.parse_args()
    verify_gate_owner()
    profile = json.loads(args.profile.read_text())
    args.out.mkdir(parents=True, exist_ok=False)
    row = dict(status='validating_contraction_runtime', started=time.time(), pid=os.getpid(),
        CPU_affinity=sorted(os.sched_getaffinity(0)), profile=str(args.profile),
        profile_sha256=digest(args.profile), program=str(args.program.resolve()))
    before = sys.argv
    try:
        row.update(validate(profile)); row['status'] = 'declared_runtime_and_FP64_smoke_passed'
        atomic_write_json(args.out / 'status.json', row)
        sys.argv = [str(args.program), *args.arguments]
        runpy.run_path(str(args.program), run_name='__main__')
        row['status'] = 'task_returned_under_checked_cutensor_runtime'
    except BaseException as error:
        row.update(status='runtime_or_task_failed_preserved', error=repr(error)); raise
    finally:
        sys.argv = before; row['finished'] = time.time()
        atomic_write_json(args.out / 'status.json', row)


if __name__ == '__main__':
    main()
