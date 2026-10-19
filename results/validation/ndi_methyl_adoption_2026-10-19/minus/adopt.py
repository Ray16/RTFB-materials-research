"""Atomically install the reviewed state; archive and verify every old byte."""
import fcntl
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import time

for key in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[key]='1'
from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity
with allocation_guard():host_capacity(1);os.sched_setaffinity(0,cpu_allocation(1))
from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.core.protocol import load_record
from redox.qm.native_thermal import load_native_thermal

BASE=Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    prepared=json.loads((BASE/'preparation.json').read_text())
    review=json.loads((BASE/'independent_review.json').read_text())
    assert review['status']=='independent_review_passed' and review['failed_checks']==0
    for path,digest in review['evidence_sha256'].items():assert sha(path)==digest,path
    owner=prepared['producer_controller']
    remote='from pathlib import Path;import sys;p=Path("/proc")/sys.argv[1]/"stat";print("live" if p.exists() and p.read_text().rsplit(")",1)[1].split()[19]==sys.argv[2] else "gone")'
    result=subprocess.run(['ssh','-o','BatchMode=yes',owner['host'],shlex.join(['python3','-c',remote,str(owner['controller_pid']),owner['process_start']])],capture_output=True,text=True,timeout=15)
    assert result.returncode==0 and result.stdout.strip()=='gone','Derivative controller must be terminal'
    case=Path(prepared['qualified_case'])
    worker=json.loads((case/'status.json').read_text())
    result=subprocess.run(['ssh','-o','BatchMode=yes',worker['host'],shlex.join(['python3','-c',remote,str(worker['pid']),worker['process_start']])],capture_output=True,text=True,timeout=15)
    assert result.returncode==0 and result.stdout.strip()=='gone','Numerical producer must be terminal'
    canonical=Path(prepared['canonical']);source=Path(prepared['prepared']);install=BASE/'install_state'
    shutil.copytree(source,install)
    spec=importlib.util.spec_from_file_location('promoted_state_independent',ROOT/'scripts/validation/audit/recompute_axes.py')
    audit=importlib.util.module_from_spec(spec);spec.loader.exec_module(audit)
    q,m=prepared['charge'],prepared['mult']
    def verify(directory):
        native=load_native_thermal(directory,q,m);energy=load_record(directory,q,m)
        raw=json.loads((directory/'result.json').read_text())
        correction,imaginary,definition=audit.thermal_inputs(directory,raw)
        assert native and energy and correction is not None and imaginary==0
        assert definition=='native_solution_dft_qrrho_spin' and abs(correction-native['g_thermal_total_eV'])<1e-8
        return correction
    verify(install)
    locks=ROOT/'logs/pipeline_locks'
    with (locks/(prepared['id']+'.lock')).open('a+') as general,(locks/(prepared['id']+'_native_thermal.lock')).open('a+') as native_lock:
        for lock in (general,native_lock):fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        for path,digest in prepared['canonical_inputs_sha256'].items():assert sha(path)==digest,path
        hashes={str(p.relative_to(canonical)):sha(p) for p in canonical.rglob('*') if p.is_file()}
        atomic_write_json(BASE/'superseded_SHA256.json',hashes)
        archive=BASE/'superseded_canonical'
        assert canonical.stat().st_dev==BASE.stat().st_dev and not archive.exists()
        os.replace(canonical,archive)
        try:
            os.replace(install,canonical)
            correction=verify(canonical)
            assert all(sha(archive/path)==digest for path,digest in hashes.items())
        except BaseException:
            if canonical.exists():os.replace(canonical,BASE/'failed_install_preserved')
            os.replace(archive,canonical)
            raise
        receipt=dict(status='independently_reviewed_native_minimum_promoted',id=prepared['id'],state=prepared['state'],
                     geometry_sha256=sha(canonical/'opt.xyz'),old_state_archive=str(archive),
                     old_bytes_verified=True,independent_checks=review['checks'],
                     minimum_frequency_cm=review['minimum_frequency_cm'],corrected_thermal_eV=correction,
                     source_review=str(BASE/'independent_review.json'),reviewed_state_archive=str(source),
                     source_producer_controller_terminal=True,source_numerical_worker_terminal=True,
                     dependent_cross_energies='Retained as addressed history; changed geometry invalidates their current use.',
                     experimental_accuracy_validated=False,finished=time.time())
        atomic_write_json(BASE/'adoption_receipt.json',receipt)
        print(json.dumps(receipt,indent=2))


if __name__=='__main__':main()
