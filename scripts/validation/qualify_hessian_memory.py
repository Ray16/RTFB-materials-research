#!/usr/bin/env python
"""Replay a certified singlet/radical Hessian with a larger execution budget."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import sys
import time
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'src'))
from redox.core.common import DFT
from redox.qm.native_thermal import ensure_native_thermal,load_native_thermal,cache_path,NATIVE_THERMAL
from redox.qm.hessian_memory import certificate,ALLOCATION_POLICY
from redox.validation.spins import atomic_json


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def run_control(out,state):
    gid='exp_naphthalimide_control';sd=DFT/gid/state
    raw=json.loads((sd/'result.json').read_text());q,m=raw['charge'],raw['mult']
    baseline=load_native_thermal(sd,q,m)
    if baseline is None:raise RuntimeError('Qualification requires a valid six-GiB baseline')
    baseline_path,_=cache_path(sd,q,m)
    baseline_h=baseline_path.with_suffix('')/'hessian_cartesian.npy'
    work=out/state;work.mkdir()
    for name in ['opt.xyz','result.json']:shutil.copy2(sd/name,work/name)
    print('MEMORY_CONTROL_START',state,flush=True)
    _,trial=ensure_native_thermal(gid,state,state_dir=work,memory_cap_gib=24,qualification=True)
    trial_path,_=cache_path(work,q,m);trial_h=trial_path.with_suffix('')/'hessian_cartesian.npy'
    A=np.load(baseline_h);B=np.load(trial_h)
    rel=float(np.linalg.norm(A-B)/np.linalg.norm(A));maximum=float(abs(A-B).max())
    thermal=max(abs(baseline['nuclear_by_cutoff_eV'][str(k)]-trial['nuclear_by_cutoff_eV'][str(k)]) for k in (25,50,100))
    passed=rel<=1e-6 and maximum<=1e-6 and thermal<=1e-6 and sha(sd/'opt.xyz')==sha(work/'opt.xyz')
    evidence=[baseline_path,baseline_h,sd/'opt.xyz',trial_path,trial_h,work/'opt.xyz']
    report=dict(status='control_verified' if passed else 'control_failed',state=state,mult=m,
                relative_matrix_error=rel,maximum_matrix_difference=maximum,
                maximum_thermal_difference_eV=thermal,native_gates_passed=passed,
                evidence_sha256={str(p.relative_to(ROOT)):sha(p) for p in evidence},finished=time.time())
    atomic_json(out/(state+'_control.json'),report)
    print('MEMORY_CONTROL_RESULT',json.dumps(report),flush=True)
    if not passed:raise RuntimeError('Execution budget replay differs from the qualified baseline')


def seal(out):
    controls=[json.loads((out/(s+'_control.json')).read_text()) for s in ('neu','red1')]
    if any(c['status']!='control_verified' for c in controls):raise RuntimeError('Both controls must pass')
    data=dict(status='execution_budget_qualified',schema=1,protocol=NATIVE_THERMAL,
              software={'pyscf':'2.14.0','gpu4pyscf':'1.8.1'},caps_GiB=[24],
              scientific_model_changed=False,controls=controls,created=time.time(),allocation_policy=ALLOCATION_POLICY)
    data['content_sha256']=hashlib.sha256(json.dumps(data,sort_keys=True,separators=(',',':')).encode()).hexdigest()
    path=ROOT/'config/native_hessian_memory_validation.json'
    if path.exists():raise RuntimeError('Existing memory certificate must be preserved and reviewed')
    atomic_json(path,data)
    assert certificate(NATIVE_THERMAL)==data['content_sha256']
    atomic_json(out/'qualification_status.json',dict(status='execution_budget_qualified',certificate=str(path),sha256=sha(path)))
    print('MEMORY_BUDGET_QUALIFIED',flush=True)


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True)
    ap.add_argument('--state',choices=['neu','red1']);ap.add_argument('--seal',action='store_true');a=ap.parse_args()
    a.out.mkdir(parents=True,exist_ok=True)
    if a.seal:seal(a.out)
    else:
        if not a.state or not os.getenv('CUDA_VISIBLE_DEVICES'):ap.error('control requires a reserved GPU and state')
        run_control(a.out,a.state)


if __name__=='__main__':main()
