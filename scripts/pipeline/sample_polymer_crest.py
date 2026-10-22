#!/usr/bin/env python
"""Run a durable CPU CREST search and retain exact-graph/stereochemistry checks."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import signal
import subprocess
import time

from redox.core.io import atomic_write_json
from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity


def xyz_frames(path):
    with Path(path).open() as stream:
        while True:
            first=stream.readline()
            if not first:return
            n=int(first.strip());comment=stream.readline()
            lines=[stream.readline() for _ in range(n)]
            if n<1 or not comment or any(not line for line in lines):
                raise ValueError('Truncated CREST geometry must not enter the conformer pool')
            yield str(n)+'\n'+comment+''.join(lines)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input',type=Path,required=True);parser.add_argument('--case',type=Path,required=True)
    parser.add_argument('--out',type=Path,required=True);parser.add_argument('--dry',action='store_true')
    parser.add_argument('--threads',type=int,default=2);parser.add_argument('--max-hours',type=float,default=24)
    parser.add_argument('--crest',type=Path,default=Path('/nfs/lambda_stor_01/homes/rzhu/miniforge3/envs/crest/bin/crest'))
    parser.add_argument('--energy-window-kcal',type=float,default=8.)
    args=parser.parse_args()
    if not 0<args.max_hours<=48 or args.threads<1 or args.energy_window_kcal<=0:
        parser.error('Positive resource settings and a finite window <=48 hours required')
    with allocation_guard():
        capacity=host_capacity(args.threads);os.sched_setaffinity(0,cpu_allocation(args.threads))
    caps=('OMP_NUM_THREADS','OMP_THREAD_LIMIT','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS','XTB_THREADS')
    if any(os.environ.get(k)!=str(args.threads) for k in caps):
        raise ValueError('Native thread caps must match the assigned CPU slots')
    args.out.mkdir(parents=True,exist_ok=False)
    case=json.loads(args.case.read_text())
    from redox.build.geometry import geometry_stereo_status
    from redox.validation.protonation import hydrogen_topology_status
    if hydrogen_topology_status(case['smiles'],args.input)!='ok' or geometry_stereo_status(case['smiles'],args.input) not in ('ok','not_applicable'):
        raise ValueError('Input is not the exact declared stereochemical graph')
    seed=args.out/'input.xyz';seed.write_bytes(args.input.read_bytes())
    xtb=args.crest.with_name('xtb')
    argv=[str(args.crest.resolve()),'input.xyz','--gfn2','--alpb','acetonitrile',
          '--T',str(args.threads),'--ewin',str(args.energy_window_kcal),
          '--keepdir','--nosz','-xnam',str(xtb.resolve())]
    if args.dry:argv.append('--dry')
    sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
    source={str(p.resolve()):sha(p) for p in (Path(__file__),args.input,args.case,args.crest,xtb)}
    status=dict(status='starting',host=os.uname().nodename,pid=os.getpid(),
        process_start=Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()[19],
        CPU_affinity=sorted(os.sched_getaffinity(0)),thread_caps={k:os.environ[k] for k in caps},
        capacity=capacity,argv=argv,source_sha256=source,case=case,started=time.time(),
        low_level_sampling_only=True,DFT_reranking_complete=False,ensemble_G_computed=False,
        MD_seed='Pinned CREST 3.0.2 program default; no documented iMTD seed override. Realized trajectories retained.')
    atomic_write_json(args.out/'status.json',status)
    child=None
    def stop(sig,frame):raise KeyboardInterrupt(str(sig))
    signal.signal(signal.SIGTERM,stop)
    env=dict(os.environ,PATH=str(args.crest.parent)+os.pathsep+os.environ.get('PATH',''))
    try:
        with (args.out/'crest.log').open('w') as log:
            child=subprocess.Popen(argv,cwd=args.out,env=env,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
            status.update(status='running',predictor_pid=child.pid,
                predictor_start=Path(f'/proc/{child.pid}/stat').read_text().rsplit(')',1)[1].split()[19])
            atomic_write_json(args.out/'status.json',status)
            rc=child.wait(timeout=args.max_hours*3600)
        if rc:raise RuntimeError('CREST did not finish successfully; all outputs preserved')
        if any(sha(p)!=h for p,h in source.items()):raise ValueError('Sampling source/input changed during execution')
        if args.dry:
            status.update(status='dry_settings_check_complete_no_sampling')
        else:
            ensemble=args.out/'crest_conformers.xyz'
            if not ensemble.exists():raise ValueError('Successful process produced no conformer ensemble')
            output=args.out/'checked_conformers';output.mkdir();records=[]
            for index,frame in enumerate(xyz_frames(ensemble)):
                path=output/f'conf_{index:05d}.xyz';path.write_text(frame)
                graph=hydrogen_topology_status(case['smiles'],path)
                stereo=geometry_stereo_status(case['smiles'],path)
                records.append(dict(index=index,path=str(path),sha256=sha(path),graph_QC=graph,stereo_QC=stereo,
                    eligible_for_DFT_reranking=graph=='ok' and stereo in ('ok','not_applicable')))
            atomic_write_json(args.out/'conformers.json',dict(records=records,ensemble_sha256=sha(ensemble)))
            status.update(status='low_level_ensemble_complete_DFT_pending',conformers=len(records),
                topology_and_stereo_pass=sum(r['eligible_for_DFT_reranking'] for r in records))
    except BaseException as error:
        status.update(status='failed_or_interrupted_preserving_outputs',error=repr(error));raise
    finally:
        if child is not None and child.poll() is None:
            os.killpg(child.pid,signal.SIGTERM)
            try:child.wait(timeout=30)
            except subprocess.TimeoutExpired:os.killpg(child.pid,signal.SIGKILL);child.wait()
        status['finished']=time.time();atomic_write_json(args.out/'status.json',status)


if __name__=='__main__':main()
