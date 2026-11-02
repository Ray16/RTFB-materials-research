#!/usr/bin/env python3
"""Repeat the unchanged primary native implementation under corrected screening."""
import argparse
import json
import math
import os
from pathlib import Path
import time
from importlib.util import spec_from_file_location,module_from_spec


def load(path,name):
    s=spec_from_file_location(name,path);m=module_from_spec(s);s.loader.exec_module(m);return m


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);a=p.parse_args()
    from redox.core.common import ROOT
    from redox.core.io import atomic_write_json
    from redox.workflows.resources import verify_gate_owner
    verify_gate_owner();case=a.plan.parent;sources=case.parent/'sources';plan=json.loads(a.plan.read_text())
    sha256_file=load(sources/'hashing.py','primary_native_frozen_hashing').sha256_file
    if (plan.get('role')!='primary_molecular_Hessian_control' or plan.get('geometry_repair_permitted')
            or any(sha256_file(p)!=h for p,h in plan['source_sha256'].items())):
        raise ValueError('Source-exact fixed-geometry primary method control required')
    if (len(os.sched_getaffinity(0))!=2 or any(os.environ.get(k)!='2' for k in
            ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS','XTB_THREADS'))):
        raise ValueError('Two disjoint CPUs and matching six native caps required')
    owner=dict(host=os.uname().nodename,pid=os.getpid(),process_start=Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()[19])
    record=dict(status='primary_molecular_Hessian_control_running',owner=owner,started=time.time(),native_G_qualified=False)
    atomic_write_json(case/'status.json',record)
    try:
        from redox.qm import dft
        from redox.qm.native_thermal import ensure_native_thermal,cache_path
        helper=load(sources/'helper.py','primary_native_screening_helper')
        checkpoint=load(sources/'checkpoints.py','primary_native_checkpoints')
        monitor=load(ROOT/'scripts/pipeline/run_pi_stacking.py','primary_native_GPU_monitor')
        dft._build_smd_mf=checkpoint.recording_builder(dft._build_smd_mf,case/'SCF_evaluations',owner)
        directory=case/'state'
        with helper.corrected_screening(sources/'kernel.cu',case/'screening_implementation'), monitor.MemoryMonitor() as memory:
            kind,thermal=ensure_native_thermal(plan['native_id'],plan['native_state'],state_dir=directory,
                memory_cap_gib=plan['native_memory_cap_GiB'])
        if kind!='computed' or not 0<memory.peak<=plan['maximum_peak_GPU_MiB']:
            raise ValueError('Fresh native execution/resource check failed')
        native_path,_=cache_path(directory,plan['charge'],plan['mult'])
        energy=json.loads((sources/'corrected_electronic.json').read_text())
        values={c:energy['e_smd_eV']+thermal['nuclear_by_cutoff_eV'][c]+thermal['g_spin_partition_eV'] for c in ('25','50','100')}
        if not all(math.isfinite(v) for v in values.values()):raise ValueError('Finite composite component values required')
        atomic_write_json(case/'result.json',dict(key=plan['source_id'],charge=plan['charge'],mult=plan['mult'],
            corrected_electronic_eV=energy['e_smd_eV'],native_record_path=str(native_path),
            nuclear_by_cutoff_eV=thermal['nuclear_by_cutoff_eV'],spin_partition_eV=thermal['g_spin_partition_eV'],
            G_by_cutoff_eV=values,reference_G_difference_eV={c:values[c]-plan['reference_checked_G_by_cutoff_eV'][c] for c in values},
            geometry_repaired=False,native_G_qualified=False,material_admission=False))
        if any(sha256_file(p)!=h for p,h in plan['source_sha256'].items()):raise ValueError('Native control source changed')
        record.update(status='primary_molecular_Hessian_control_finished_review_pending',observed_peak_GPU_MiB=memory.peak)
    except BaseException as e:
        record.update(status='primary_molecular_Hessian_control_failed_preserved',error=repr(e));raise
    finally:
        record['finished']=time.time();atomic_write_json(case/'status.json',record)


if __name__=='__main__':main()
