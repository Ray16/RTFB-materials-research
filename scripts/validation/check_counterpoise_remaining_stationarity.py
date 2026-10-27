#!/usr/bin/env python3
"""Independently check one reviewed remaining CPU counterpoise term."""
import argparse
import json
import os
from pathlib import Path
import signal
import threading
import time

import numpy as np
import psutil
import pyscf
from pyscf import lib

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api
from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity

EVALUATOR=ROOT/'scripts/validation/check_counterpoise_remaining_terms.py'
CONTROL=ROOT/'results/reviews/2026-10-26/counterpoise_cpu_pilot/remaining_water_check.json'
EXPECTED={'AB','B_own','B_ghost'}


def check(raw,review_path,out):
    raw=raw.resolve();review_path=review_path.resolve();out=out.resolve()
    if out.exists():raise FileExistsError('Preserve prior stationarity attempts')
    if not raw.is_relative_to((ROOT/'calcs').resolve()) or not out.is_relative_to((ROOT/'calcs').resolve()):raise ValueError('Use durable project storage')
    scratch=Path(os.environ.get('PYSCF_TMPDIR','/nonexistent')).resolve()
    if not scratch.is_dir() or not scratch.is_relative_to((ROOT/'calcs').resolve()):raise ValueError('Durable PySCF scratch required')
    review=json.loads(review_path.read_text())
    if (review['status']!='CPU_remaining_term_checkpoint_and_sampled_resources_reviewed_gradient_pending'
            or review['physical_gradient_recomputed'] or review['full_counterpoise_computed']):
        raise ValueError('Completed independent term review required')
    api.verify_sources(review['source_sha256']);control=json.loads(CONTROL.read_text());api.verify_sources(control['source_sha256'])
    if (control['status']!='remaining_water_counterpoise_terms_stationary' or not control['passed']
            or len(control['terms'])!=3 or {x['role'] for x in control['terms']}!=EXPECTED or not all(x['passed'] for x in control['terms'])):
        raise ValueError('Source-bound remaining water controls required')
    launch=json.loads((raw/'launch.json').read_text());role=launch['role']
    if role not in EXPECTED or review['role']!=role:raise ValueError('Review is for a different term')
    record=json.loads((raw/role/'energy.json').read_text());api.verify_sources(record['source_sha256'])
    if (launch['versions']!={'pyscf':pyscf.__version__,'numpy':np.__version__}
            or record['status']!='CPU_remaining_term_converged_independent_review_pending'
            or record['role']!=role or abs(record['E_Ha']-review['energy_Ha'])>1e-9):
        raise ValueError('Term numerical environment or reviewed energy differs')
    reviewer=api.load_module('remaining_stationarity_owner',ROOT/'scripts/validation/review_counterpoise_cpu_pilot.py')
    reviewer.verify_owner_stopped(record['owner'])
    plan=api.read_plan(Path(launch['case']));geometry=json.loads(Path(plan['input_geometry']).read_text())
    checkpoint=raw/role/'SCF.h5';coords=raw/'grid_coordinates_AU.npy';weights=raw/'grid_weights.npy'
    if (api.digest(checkpoint)!=record['checkpoint_sha256'] or api.digest(coords)!=record['common_grid_hashes']['grid_coordinates_AU.npy']
            or api.digest(weights)!=record['common_grid_hashes']['grid_weights.npy']):raise ValueError('Term checkpoint or common grid changed')
    bound=[raw/'launch.json',raw/'status.json',raw/role/'energy.json',checkpoint,coords,weights]
    if any(review['source_sha256'].get(str(p))!=api.digest(p) for p in bound):raise ValueError('Review does not bind this exact term checkpoint')
    paths=bound+[review_path,Path(launch['case'])/'plan.json',Path(plan['input_geometry']),EVALUATOR,CONTROL,Path(__file__).resolve()]
    sources={str(p):api.digest(p) for p in paths};out.mkdir(parents=True);started=time.time()
    process=psutil.Process();peak=[process.memory_info().rss];stop=threading.Event();limit=32*1024**3
    def monitor():
        while not stop.wait(2):
            rss=process.memory_info().rss;peak[0]=max(peak[0],rss)
            if rss>limit:os.kill(os.getpid(),signal.SIGTERM);return
    def interrupted(signum,frame):raise MemoryError('Term stationarity stopped by sampled RSS guard or termination')
    signal.signal(signal.SIGTERM,interrupted);watcher=threading.Thread(target=monitor,daemon=True);watcher.start()
    try:
        evaluator=api.load_module('remaining_term_physical_evaluator',EVALUATOR)
        computed=evaluator.evaluate(geometry,role,checkpoint,coords,weights,record['E_Ha'])
    except BaseException as error:
        atomic_write_json(out/'failure.json',dict(status='remaining_stationarity_failed_preserved',error=repr(error),
            sampled_peak_RSS_bytes=peak[0],source_sha256=sources,qualified=False));raise
    finally:stop.set();watcher.join(timeout=3)
    api.verify_sources(sources);api.verify_sources(record['source_sha256']);api.verify_sources(review['source_sha256']);api.verify_sources(control['source_sha256'])
    passed=(computed['energy_error_Ha']<=1e-8 and computed['checkpoint_energy_error_Ha']<=1e-8
            and computed['orbital_gradient_norm_au']<=computed['orbital_gradient_threshold_au'])
    result=dict(status='CPU_remaining_term_energy_and_orbital_stationarity_checked' if passed else 'remaining_stationarity_failed_preserved',
        passed=passed,**computed,sampled_peak_RSS_bytes=peak[0],seconds=time.time()-started,source_sha256=sources,
        full_counterpoise_computed=False,solvent_counterpoise_computed=False,
        scope='One frozen electronic term; no geometry gradient or complete interaction')
    atomic_write_json(out/'summary.json',result)
    if not passed:raise ValueError('Term checkpoint fails the declared stationarity threshold')
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('raw','review','out'):p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--threads',type=int,default=8);a=p.parse_args()
    if a.threads<1 or pyscf.__version__!='2.14.0':raise ValueError('Qualified PySCF and positive thread count required')
    with allocation_guard():host_capacity(a.threads,minimum_ram_gib=80);os.sched_setaffinity(0,cpu_allocation(a.threads))
    if any(os.environ.get(k)!=str(a.threads) for k in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','XTB_THREADS')):
        raise ValueError('Declared CPU threads and capped native pools required')
    lib.num_threads(a.threads);result=check(a.raw,a.review,a.out);print(json.dumps(result))
