#!/usr/bin/env python3
"""Install the explicit xTB runtime in a fresh temporary prefix and test it.

Uses the existing Python environment; this qualifies the relocated xTB runtime,
not a complete Python/CUDA deployment. Original xTB paths can be masked with
bubblewrap during the real thermodynamic checks.
"""
import argparse
import json
import os
import re
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

from redox.core.common import ROOT
from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json
from redox.qm import qrrho
from redox.qm.xtb_runtime import validate_prefix


PROGRAM='''import json,sys,subprocess,numpy as np
from redox.qm import qrrho as q,native_thermal,projected_thermal
from pathlib import Path
prefix=Path(sys.argv[1]);assert q.LIBRARY==prefix/'lib/libxtb.so.6.7.1'
assert native_thermal.XTB==projected_thermal.XTB==prefix/'bin/xtb'
moments=np.array([15649.5560990126,484991.786527778,485024.382510520])
frequencies=[117.29518854,161.91929949,164.05816261,218.07506352,402.96868229,558.86387616,1523.62725135,1561.31279119,3460.43559680,3633.65239963,3636.81476002,3665.34271369]
rot=.5*q.AUTORCM/moments[::-1];average=moments.mean()/((1/9.10938356e-31)*(1/.52917726)**2*1e20)
values=q.kernel_thermal(frequencies,rot,average,2*15.99940492+4*1.00794075,1,cutoff_cm=100)
expected=dict(G_Ha=.017417501635698,H_Ha=.050275771916811,TS_Ha=.032858270281112)
assert all(abs(values[k]-v)<1e-12 for k,v in expected.items())
masses=q.atomic_masses([1,1]);stretch=np.array([-1.,0.,0.,1.,0.,0.])/np.sqrt(2)
matrix=(.391/q.AUTORCM)**2*masses[0]*q.AMUTOAU*np.outer(stretch,stretch)
soft=q.nuclear_thermals(matrix,[[-.5,0.,0.],[.5,0.,0.]],[1,1],2)
assert soft['rigid_rank']==5 and soft['internal_mode_count']==1 and abs(soft['frequencies_cm'][0]-.391)<1e-10
version=subprocess.run([str(native_thermal.XTB),'--version'],capture_output=True,text=True,timeout=30,check=True)
print(json.dumps(dict(library=str(q.LIBRARY),executable=str(native_thermal.XTB),upstream_kernel=values,upstream_expected=expected,soft_frequency_cm=soft['frequencies_cm'][0],soft_mode_retained=True,kernel_execution=q.kernel_execution(),xtb_version_stdout=version.stdout,CUDA_imported=any(k.startswith(('cupy','gpu4pyscf')) for k in sys.modules))))
'''


def check(lock,manager,output,offline,bwrap=None):
    if output.exists():raise FileExistsError('Preserve the previous installation check')
    output.mkdir(parents=True);inventory=json.loads((lock/'packages.json').read_text())
    temp=Path(tempfile.mkdtemp(prefix='redox_xtb_installation_'));prefix=temp/'environment'
    result=dict(status='fresh_xtb_installation_check_started',offline=offline,original_xtb_files_masked=bwrap is not None,
        lock_sha256={str(p):sha256_file(p) for p in [lock/'packages.json',lock/'environment.explicit.txt']},
        producer_sha256=sha256_file(__file__),GPU_job_started=False,full_quantum_environment_qualified=False)
    try:
        if shutil.disk_usage(temp).free<2*2**30:raise RuntimeError('Two GiB temporary-disk headroom required')
        env=dict(os.environ,CONDA_NO_PLUGINS='true')
        argv=[str(manager),'create','--yes','--solver','classic','--no-default-packages','--prefix',str(prefix),'--file',str(lock/'environment.explicit.txt')]
        if offline:argv.append('--offline')
        with (output/'install.log').open('x') as stream:
            run=subprocess.run(argv,env=env,stdout=stream,stderr=subprocess.STDOUT,timeout=600)
        if run.returncode:raise RuntimeError('Explicit xTB installation failed; see install.log')
        installed={json.loads(p.read_text())['name']:json.loads(p.read_text()) for p in (prefix/'conda-meta').glob('*.json')}
        wanted={row['name']:row for row in inventory['packages']}
        if set(installed)!=set(wanted) or any(any(installed[name].get(k)!=row[k] for k in ('version','build','sha256')) for name,row in wanted.items()):
            raise ValueError('Installed conda artifact metadata differs from the explicit closure')
        qualified=validate_prefix(prefix,qrrho.LIBRARY_SHA256)
        atomic_write_json(output/'installed_packages.json',{name:{k:row.get(k) for k in ('version','build','url','sha256','md5')} for name,row in installed.items()})
        env.update(REDOX_XTB_PREFIX=str(prefix),LD_LIBRARY_PATH=str(prefix/'lib'),PYTHONPATH=str(ROOT/'src'),TMPDIR=str(temp),PYTHONDONTWRITEBYTECODE='1')
        env.pop('REDOX_XTB_LIBRARY',None);env.pop('REDOX_XTB_EXECUTABLE',None)
        linkage={}
        for name in ('library','binary'):
            run=subprocess.run(['ldd',qualified[name]],env=env,capture_output=True,text=True,timeout=30,check=True)
            if 'not found' in run.stdout:raise ValueError('Fresh xTB prefix has an unresolved dynamic dependency')
            paths=re.findall(r'=>\s+(/\S+)',run.stdout)
            if any(not any(Path(path).resolve().is_relative_to(root) for root in (prefix.resolve(),Path('/lib').resolve(),Path('/lib64').resolve(),Path('/usr/lib').resolve())) for path in paths):
                raise ValueError('Fresh xTB linkage borrowed a non-system library outside its locked prefix')
            linkage[name]=run.stdout
        command=[sys.executable,'-c',PROGRAM,str(prefix)]
        if bwrap is not None:
            original=inventory['runtime']
            command=[str(bwrap),'--ro-bind','/','/','--bind',str(temp),str(temp),'--ro-bind',str(prefix),str(prefix),
                '--ro-bind','/dev/null',original['library'],'--ro-bind','/dev/null',original['binary'],
                '--proc','/proc','--dev-bind','/dev','/dev','--chdir',str(ROOT),'--',*command]
        run=subprocess.run(command,env=env,capture_output=True,text=True,timeout=120)
        (output/'thermo.stdout').write_text(run.stdout);(output/'thermo.stderr').write_text(run.stderr)
        if run.returncode:raise RuntimeError('Relocated thermodynamic fixture failed; see thermo.stderr')
        numeric=json.loads(run.stdout)
        if numeric['CUDA_imported']:raise ValueError('The CPU runtime check unexpectedly imported CUDA modules')
        result.update(status='fresh_locked_xtb_runtime_and_relocated_thermo_checked',installed_packages=len(installed),
            qualified_runtime=qualified,dynamic_linkage=linkage,numerical_checks=numeric,Python_environment_recreated=False)
    except BaseException as error:
        result.update(status='fresh_xtb_installation_check_failed_preserved',error=repr(error));raise
    finally:
        atomic_write_json(output/'summary.json',result);shutil.rmtree(temp)
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--lock',type=Path,required=True);p.add_argument('--manager',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True);p.add_argument('--offline',action='store_true');p.add_argument('--bwrap',type=Path);a=p.parse_args()
    from redox.workflows.resources import allocation_guard,host_capacity,cpu_allocation
    with allocation_guard():host_capacity(1,minimum_ram_gib=4);os.sched_setaffinity(0,cpu_allocation(1))
    for name in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS','XTB_THREADS'):os.environ[name]='1'
    result=check(a.lock.resolve(),a.manager.resolve(),a.out.resolve(),a.offline,None if a.bwrap is None else a.bwrap.resolve())
    print(json.dumps({k:result[k] for k in ('status','installed_packages','original_xtb_files_masked')}))
