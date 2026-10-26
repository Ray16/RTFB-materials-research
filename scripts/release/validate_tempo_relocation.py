#!/usr/bin/env python3
"""Extract the TEMPO addendum and replay it with original-project reads blocked."""
import argparse
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import runpy
import shutil
import subprocess
import sys
import tarfile
import time


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def worker(package,out,forbidden,full_native,expected_prefix):
    blocked=[]
    def audit(event,args):
        if event=='open' and args and isinstance(args[0],(str,bytes,os.PathLike)):
            path=Path(os.fsdecode(args[0])).resolve()
            if path.is_relative_to(forbidden):
                blocked.append(str(path));raise RuntimeError('Original project reads are forbidden')
    sys.addaudithook(audit)
    target=out/('native.json' if full_native else 'stdlib.json')
    sys.argv=[str(package/'reproduce.py'),'--data',str(package),'--out',str(target)]
    if full_native:sys.argv.append('--full-native')
    runpy.run_path(str(package/'reproduce.py'),run_name='__main__')
    forbidden_modules={'redox','pyscf','gpu4pyscf','cupy','torch','ase'}
    if not full_native:forbidden_modules.add('numpy')
    loaded=[n for n in sys.modules if n.split('.')[0] in forbidden_modules]
    if blocked or loaded:raise RuntimeError('Replay attempted original reads or loaded forbidden modules')
    numpy=sys.modules.get('numpy')
    if full_native and (numpy is None or numpy.__version__!='2.2.6' or Path(sys.prefix).resolve()!=expected_prefix
            or not Path(numpy.__file__).resolve().is_relative_to(expected_prefix)):
        raise ValueError('Full replay did not use the locked NumPy environment')
    result=json.loads(target.read_text())
    if result['status']!='portable_TEMPO_local_oxidation_reconstructed' or result['states']!=4:
        raise ValueError('Unexpected TEMPO replay result')
    if full_native and (result['native_reconstruction']['Hessians_reconstructed']!=4 or not result['nuclear_partition_reintegrated']):
        raise ValueError('Four independent Hessian reconstructions required')
    proof=dict(status='TEMPO_relocated_full_native_checked' if full_native else 'TEMPO_relocated_stdlib_checked',
        checked=time.time(),original_project_open_attempts=blocked,forbidden_modules_loaded=loaded,
        interpreter=sys.executable,isolated_python=bool(sys.flags.isolated),numpy_version=numpy.__version__ if numpy else None,
        result_sha256=sha(target),validator_sha256=sha(__file__))
    (out/('native_isolation.json' if full_native else 'stdlib_isolation.json')).write_text(json.dumps(proof,indent=2)+'\n')


def validate(source,archive,out,forbidden,stdlib_python,fresh_python,expected_prefix):
    if out.exists():raise FileExistsError('Preserve previous relocation attempts')
    out.mkdir(parents=True)
    tar_path=out/'payload.tar'
    with tar_path.open('xb') as stream:
        subprocess.run(['zstd','-dc',str(archive)],stdout=stream,check=True,timeout=120)
    expected_root=source.name
    with tarfile.open(tar_path,'r:') as tf:
        members=tf.getmembers()
        for member in members:
            path=PurePosixPath(member.name)
            if (path.is_absolute() or '..' in path.parts or not path.parts or path.parts[0]!=expected_root
                    or not (member.isdir() or member.isfile())):
                raise ValueError('Unsafe archive member')
        tf.extractall(out,filter='data')
    package=out/expected_root
    if sha(package/'manifest.json')!=sha(source/'manifest.json'):
        raise ValueError('Extracted manifest differs from staged source')
    manifest=json.loads((package/'manifest.json').read_text())
    expected=set(manifest['files_sha256'])|{'manifest.json'}
    actual={str(p.relative_to(package)) for p in package.rglob('*') if p.is_file()}
    if actual!=expected or any(sha(package/p)!=h for p,h in manifest['files_sha256'].items()):
        raise ValueError('Extracted inventory or checksum differs')
    validator=out/'validator.py';shutil.copyfile(__file__,validator)
    env=dict(os.environ,OMP_NUM_THREADS='1',OMP_THREAD_LIMIT='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',NUMEXPR_NUM_THREADS='1')
    env.pop('LD_LIBRARY_PATH',None);env.pop('PYTHONPATH',None)
    for mode,python in [('stdlib',stdlib_python),('native',fresh_python)]:
        args=[str(python),'-I']+(['-S'] if mode=='stdlib' else [])+[str(validator),'--worker','--package',str(package),
            '--out',str(out),'--forbid-root',str(forbidden)]
        if mode=='native':args+=['--full-native','--expected-prefix',str(expected_prefix)]
        with (out/(mode+'_worker.log')).open('x') as log:
            subprocess.run(args,cwd=out,env=env,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=120)
    stdlib=json.loads((out/'stdlib.json').read_text());native=json.loads((out/'native.json').read_text())
    proof=dict(status='TEMPO_archive_and_two_relocated_replays_checked',checked=time.time(),archive=str(archive),
        archive_sha256=sha(archive),archive_bytes=archive.stat().st_size,package_manifest_sha256=sha(package/'manifest.json'),
        archived_regular_files=len(expected),states=4,Hessians_reconstructed=4,stdlib_result=stdlib,
        native_result=native,stdlib_isolation=json.loads((out/'stdlib_isolation.json').read_text()),
        native_isolation=json.loads((out/'native_isolation.json').read_text()),validator_sha256=sha(__file__),
        quantum_calculations_repeated=False,public_archive_deposited=False,experimental_accuracy_inferred=False)
    (out/'validation.json').write_text(json.dumps(proof,indent=2)+'\n')
    return proof


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--worker',action='store_true');p.add_argument('--full-native',action='store_true')
    for name in ('source','archive','package','out','forbid-root','stdlib-python','fresh-python','expected-prefix'):
        p.add_argument('--'+name,type=Path)
    a=p.parse_args()
    if a.worker:
        worker(a.package.resolve(),a.out.resolve(),a.forbid_root.resolve(),a.full_native,a.expected_prefix.resolve() if a.expected_prefix else None)
    else:
        if not all((a.source,a.archive,a.out,a.forbid_root,a.stdlib_python,a.fresh_python,a.expected_prefix)):
            p.error('Source, archive, out, forbidden root and both interpreters are required')
        result=validate(a.source.resolve(),a.archive.resolve(),a.out.resolve(),a.forbid_root.resolve(),
                        a.stdlib_python.resolve(),a.fresh_python.absolute(),a.expected_prefix.resolve())
        print(json.dumps({'status':result['status'],'states':result['states'],'Hessians':result['Hessians_reconstructed']}))
