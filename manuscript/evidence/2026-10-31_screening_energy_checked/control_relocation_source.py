#!/usr/bin/env python3
"""Verify/extract the control archive, then replay with original reads blocked."""
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


def worker(example, out, mode, forbidden, expected_prefix):
    blocked = []
    def audit(event, args):
        if event == 'open' and isinstance(args[0], (str, bytes)):
            value = args[0].decode() if isinstance(args[0], bytes) else args[0]
            path = Path(value).resolve()
            if path.is_relative_to(forbidden):
                blocked.append(str(path))
                raise RuntimeError('Opening the original project is forbidden: '+str(path))
    sys.addaudithook(audit)
    result = out / (mode+'_replay.json')
    sys.argv = [str(example/'reproduce.py'), '--package', str(example), '--out', str(result)]
    if mode == 'electronic': sys.argv.append('--electronic-only')
    runpy.run_path(str(example/'reproduce.py'), run_name='__main__')
    forbidden_modules = sorted(n for n in sys.modules
        if n.split('.')[0] in {'redox','pyscf','gpu4pyscf','cupy','torch','ase'})
    numpy = sys.modules.get('numpy')
    if blocked or forbidden_modules:
        raise RuntimeError('Replay attempted an original path or loaded a quantum/production module')
    if mode == 'full':
        if (numpy is None or numpy.__version__ != '2.2.6' or Path(sys.prefix).resolve() != expected_prefix
                or not Path(numpy.__file__).resolve().is_relative_to(expected_prefix)):
            raise RuntimeError('Full replay did not use the fresh pinned NumPy environment')
    elif numpy is not None:
        raise RuntimeError('Standard-library replay unexpectedly loaded NumPy')
    proof = dict(status='extracted_control_replay_passed_with_original_project_reads_blocked', mode=mode,
        checked=time.time(), original_project_open_attempts=blocked,
        forbidden_quantum_modules_loaded=forbidden_modules, interpreter=sys.executable,
        python_version=sys.version, CPU_affinity=sorted(os.sched_getaffinity(0)),
        numpy_version=None if numpy is None else numpy.__version__,
        numpy_location=None if numpy is None else numpy.__file__, result_sha256=sha(result),
        validator_sha256=sha(__file__), example_manifest_sha256=sha(example/'data/manifest.json'))
    (out/(mode+'_proof.json')).write_text(json.dumps(proof,indent=2)+'\n')


def validate(source, archive, fresh_python, stdlib_python, out):
    root = Path(__file__).resolve().parents[2]
    if out.is_relative_to(root):
        raise ValueError('Relocated replay must live outside the original project')
    from redox.workflows.resources import allocation_guard, cpu_allocation, host_capacity
    with allocation_guard(): host_capacity(1); os.sched_setaffinity(0,cpu_allocation(1))
    out.mkdir(parents=True,exist_ok=False)
    expected = json.loads((source/'package_manifest.json').read_text())['files_sha256']
    expected['package_manifest.json'] = sha(source/'package_manifest.json')
    extracted = out / 'extracted'; extracted.mkdir()
    members = {}; process = subprocess.Popen(['zstd','-dc',str(archive)],stdout=subprocess.PIPE)
    try:
        with tarfile.open(fileobj=process.stdout, mode='r|') as stream:
            for item in stream:
                parts = PurePosixPath(item.name).parts
                if (not parts or parts[0] != source.name or '..' in parts
                        or PurePosixPath(item.name).is_absolute() or not (item.isdir() or item.isfile())):
                    raise ValueError('Archive member is not a contained regular file/directory')
                target = extracted.joinpath(*parts)
                if item.isdir(): target.mkdir(parents=True,exist_ok=True); continue
                relative = str(PurePosixPath(*parts[1:]))
                if relative not in expected or relative in members:
                    raise ValueError('Unexpected or duplicate archived file: '+relative)
                target.parent.mkdir(parents=True,exist_ok=True)
                with stream.extractfile(item) as incoming, target.open('xb') as output:
                    shutil.copyfileobj(incoming,output)
                actual = sha(target)
                if actual != expected[relative]: raise ValueError('Archived bytes changed: '+relative)
                members[relative] = actual
        process.stdout.close()
        if process.wait() != 0 or members != expected:
            raise ValueError('Archive failed decompression or has an incomplete inventory')
    finally:
        if process.poll() is None:
            process.terminate(); process.wait()
    helper = out/'validate_source.py'; shutil.copyfile(Path(__file__).resolve(),helper)
    example = extracted/source.name/'example'
    env = dict(os.environ); env.pop('PYTHONPATH',None); env.pop('LD_LIBRARY_PATH',None)
    for key in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS'):
        env[key]='1'
    for mode, interpreter, flags in [('electronic',stdlib_python,['-I','-S']),('full',fresh_python,['-I'])]:
        argv = [str(interpreter),*flags,str(helper),'--worker',mode,'--example',str(example),
            '--out',str(out),'--forbid-root',str(root),'--expected-prefix',str(fresh_python.parent.parent)]
        with (out/(mode+'.log')).open('x') as log:
            subprocess.run(argv,env=env,cwd=out,stdout=log,stderr=subprocess.STDOUT,check=True)
    full = json.loads((out/'full_replay.json').read_text())
    if full['quantum_states_reconstructed'] != 8 or full['Fc_functional_controls']['electronic_states'] != 6:
        raise ValueError('The relocated replay did not cover the declared raw states')
    proof = dict(status='archive_inventory_and_two_extracted_replays_checked',checked=time.time(),
        archive=str(archive),archive_sha256=sha(archive),archive_bytes=archive.stat().st_size,
        archive_regular_files=len(members),package_manifest_sha256=sha(source/'package_manifest.json'),
        source_code_commit=json.loads((source/'package_manifest.json').read_text())['source_code_commit'],
        raw_quantum_control_states=8,Fc_electronic_records=6,raw_SCF_evaluations=392,baseline_evaluations=8,displaced_evaluations=384,
        isolated_replays={mode:json.loads((out/(mode+'_proof.json')).read_text()) for mode in ('electronic','full')},
        validator_sha256=sha(__file__),source_validator=str(Path(__file__).resolve()),
        public_archive_deposited=False,polymer_G_qualification_inferred=False)
    (out/'summary.json').write_text(json.dumps(proof,indent=2)+'\n')
    print(json.dumps({k:v for k,v in proof.items() if k != 'isolated_replays'}))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--worker',choices=('electronic','full'))
    for name in ('example','out','forbid-root','expected-prefix','source','archive','fresh-python','stdlib-python'):
        parser.add_argument('--'+name,type=Path)
    args=parser.parse_args()
    if args.worker:
        if any(getattr(args,k) is None for k in ('example','out','forbid_root','expected_prefix')):
            parser.error('worker requires example/out/forbid-root/expected-prefix')
        worker(args.example.resolve(),args.out.resolve(),args.worker,args.forbid_root.resolve(),args.expected_prefix.resolve())
    else:
        if any(getattr(args,k) is None for k in ('source','archive','fresh_python','stdlib_python','out')):
            parser.error('validation requires source/archive/fresh-python/stdlib-python/out')
        # Retain the declared venv entry point; resolving its symlink changes its prefix.
        validate(args.source.resolve(),args.archive.resolve(),args.fresh_python.absolute(),
                 args.stdlib_python.absolute(),args.out.resolve())


if __name__=='__main__': main()
