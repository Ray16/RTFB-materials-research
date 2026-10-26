#!/usr/bin/env python
"""Review relocated qualified xTB bindings against kernel and CLI controls."""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import time

from redox.core.common import ROOT
from redox.qm.vibrations import parse_xtb_thermal_eV
from redox.qm.xtb_runtime import QUALIFIED_XTB_BINARY_SHA256


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def compare_probes(rows):
    if len(rows)!=3:raise ValueError('Original, relocated and relocated-interpreter probes required')
    for r in rows:
        if (r['status']!='qualified_xtb_relocated_kernel_probe_checked' or not r['loaded_library_path_verified']
                or r['quantum_calculations_run'] or r['native_state_qualified'] or r['full_platform_portability_inferred']):
            raise ValueError('Probe scope or actual library verification differs')
    first=rows[0];differences=[]
    for r in rows[1:]:
        if set(r['thermal'])!=set(first['thermal']) or len(r['masses_amu'])!=len(first['masses_amu']):
            raise ValueError('Probe output shapes differ')
        differences.extend(abs(a-b) for a,b in zip(first['masses_amu'],r['masses_amu']))
        differences.extend(abs(first['thermal'][k]-r['thermal'][k]) for k in first['thermal'])
    if not all(math.isfinite(v) and v<=1e-12 for v in differences):raise ValueError('Relocated kernel output differs')
    return max(differences,default=0.)


def review(probe_paths,geometry,hessian,qualification,out):
    if out.exists():raise FileExistsError('Preserve existing relocation review')
    rows=[json.loads(p.read_text()) for p in probe_paths];difference=compare_probes(rows)
    proof=json.loads(qualification.read_text())
    if proof['status']!='fresh_compatible_GPU_deployment_independently_checked':raise ValueError('Backing control not qualified')
    for p in (geometry,hessian):
        expected=proof['original_and_fresh_backing_sha256'].get(str(p))
        if expected is None or sha(p)!=expected:raise ValueError('Geometry or Hessian does not match the checked control')
    for row in rows:
        runtime=row['runtime'];prefix=Path(runtime['prefix']).resolve()
        for label,hash_label in [('library','library_sha256'),('binary','binary_sha256')]:
            path=Path(runtime[label]).resolve()
            if not path.is_relative_to(prefix) or sha(path)!=runtime[hash_label]:raise ValueError('Runtime artifact or prefix changed')
        if (runtime['library_sha256']!=proof['libxtb_sha256'] or row['worker']['library_sha256']!=runtime['library_sha256']
                or runtime['binary_sha256']!=QUALIFIED_XTB_BINARY_SHA256):
            raise ValueError('Runtime library does not match qualification')
    if not Path(rows[2]['worker']['interpreter']).is_relative_to(Path(rows[2]['runtime']['prefix'])):
        raise ValueError('Fresh-interpreter probe did not use the restored prefix')
    out.mkdir(parents=True)
    cli=[]
    for label,row in [('original',rows[0]),('fresh',rows[2])]:
        folder=out/label;folder.mkdir();shutil.copyfile(geometry,folder/'mol.xyz');shutil.copyfile(hessian,folder/'hessian')
        env=dict(os.environ);env.pop('LD_PRELOAD',None);env['LD_LIBRARY_PATH']=str(Path(row['runtime']['library']).parent)
        for k in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS','XTB_THREADS'):env[k]='1'
        result=subprocess.run([row['runtime']['binary'],'thermo','--sthr','50','--temp','298.15','mol.xyz','hessian'],
            cwd=folder,env=env,capture_output=True,text=True,timeout=60,check=True)
        (folder/'cli.log').write_text(result.stdout+result.stderr)
        cli.append(dict(prefix=label,G_RRHO_eV=parse_xtb_thermal_eV(result.stdout+result.stderr),log_sha256=sha(folder/'cli.log')))
    cli_difference=abs(cli[0]['G_RRHO_eV']-cli[1]['G_RRHO_eV'])
    if not math.isfinite(cli_difference) or cli_difference>1e-10:raise ValueError('Relocated xTB CLI thermochemistry differs')
    source=[*probe_paths,geometry,hessian,qualification,Path(__file__).resolve(),ROOT/'src/redox/qm/xtb_runtime.py',
        ROOT/'scripts/pipeline/run_portable_thermal.py',ROOT/'src/redox/qm/qrrho.py',ROOT/'src/redox/qm/xtb_kernel.py',
        Path(rows[2]['worker']['interpreter'])]+[Path(r['runtime'][k]) for r in (rows[0],rows[2]) for k in ('binary','library')]
    result=dict(status='qualified_xtb_prefix_binding_kernel_and_CLI_relocation_checked',checked=time.time(),
        maximum_kernel_or_mass_difference=difference,CLI_difference_eV=cli_difference,CLI=cli,
        restored_prefix=rows[2]['runtime']['prefix'],restored_interpreter=rows[2]['worker']['interpreter'],
        qualified_library_sha256=rows[2]['runtime']['library_sha256'],qualified_binary_sha256=rows[2]['runtime']['binary_sha256'],
        source_sha256={str(p):sha(p) for p in source},native_GPU_entrypoint_requalified=False,
        full_platform_portability_inferred=False,quantum_calculations_run=False,original_physical_QC_repeated=False)
    (out/'summary.json').write_text(json.dumps(result,indent=2)+'\n');return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--probe',type=Path,action='append',required=True);p.add_argument('--geometry',type=Path,required=True);p.add_argument('--hessian',type=Path,required=True);p.add_argument('--qualification',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity
    with allocation_guard():host_capacity(1);os.sched_setaffinity(0,cpu_allocation(1))
    result=review([p.resolve() for p in a.probe],a.geometry.resolve(),a.hessian.resolve(),a.qualification.resolve(),a.out.resolve())
    print(json.dumps({k:result[k] for k in ('status','maximum_kernel_or_mass_difference','CLI_difference_eV')}))
