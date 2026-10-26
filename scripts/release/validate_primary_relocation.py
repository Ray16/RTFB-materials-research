#!/usr/bin/env python3
"""Reconstruct the primary package Hessians with original-project reads blocked."""
import argparse
import hashlib
from importlib import metadata
import json
import os
from pathlib import Path
import runpy
import sys
import time


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def validate_versions(found):
    expected={'numpy':'2.2.6','scipy':'1.17.1','ase':'3.29.0'}
    if found != expected:raise ValueError('Primary replay dependencies differ from the declared lock')
    return expected


def worker(package,out,forbidden,expected_prefix):
    if not sys.flags.isolated or Path(sys.prefix).resolve()!=expected_prefix.resolve():
        raise ValueError('Use the declared isolated replay environment')
    blocked=[]
    def audit(event,args):
        if event=='open' and args and isinstance(args[0],(str,bytes,os.PathLike)):
            path=Path(os.fsdecode(args[0])).resolve()
            if path.is_relative_to(forbidden):
                blocked.append(str(path));raise RuntimeError('Original project reads are forbidden')
    sys.addaudithook(audit)
    versions=validate_versions({name:metadata.version(name) for name in ('numpy','scipy','ase')})
    for name in ('numpy','scipy','ase'):
        __import__(name)
    result_path=out/'primary_hessians.json'
    program=package/'code/scripts/release/review_package_hessians.py'
    sys.argv=[str(program),str(package),'--out',str(result_path)]
    runpy.run_path(str(program),run_name='__main__')
    result=json.loads(result_path.read_text())
    if (result['status']!='packaged_Hessians_modes_and_nuclear_sums_reconstructed'
            or result['checked']!=242 or len(result['restart_geometry_frequency_checks'])!=2
            or result['quantum_calculations_run'] or result['xTB_library_loaded']
            or result['original_physical_QC_repeated'] or result['experimental_accuracy_inferred']):
        raise ValueError('Primary Hessian replay coverage or scope differs')
    forbidden_modules={'redox','pyscf','gpu4pyscf','cupy','torch'}
    loaded=[name for name in sys.modules if name.split('.')[0] in forbidden_modules]
    if blocked or loaded:raise ValueError('Replay attempted original reads or loaded quantum modules')
    for name in ('numpy','scipy','ase'):
        module=sys.modules.get(name)
        if module is None or not Path(module.__file__).resolve().is_relative_to(expected_prefix.resolve()):
            raise ValueError('Replay dependency is outside the declared environment: '+name)
    proof=dict(status='primary_242_Hessians_reconstructed_after_relocation',checked=time.time(),
        package_manifest_sha256=sha(package/'package_manifest.json'),result_sha256=sha(result_path),
        validator_sha256=sha(__file__),dependency_versions=versions,interpreter=sys.executable,
        isolated_python=True,original_project_open_attempts=blocked,forbidden_modules_loaded=loaded,
        quantum_calculations_repeated=False,original_physical_QC_repeated=False,
        experimental_accuracy_inferred=False)
    (out/'primary_isolation.json').write_text(json.dumps(proof,indent=2)+'\n')
    return proof


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('package','out','forbid-root','expected-prefix'):p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args();worker(a.package.resolve(),a.out.resolve(),a.forbid_root.resolve(),a.expected_prefix.resolve())
