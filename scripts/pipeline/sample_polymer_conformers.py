#!/usr/bin/env python
"""Checkpointed CREST sampling with bounded outer parallelism and serial BLAS.

The earlier two-core pilot producer remains unchanged for its live searches.
This successor uses the same GFN2/ALPB(MeCN), energy window and topology checks.
Only execution parallelism differs; trajectories and all numerical inputs persist.
"""
import argparse,hashlib,importlib.util,json,os,signal,subprocess,time
from pathlib import Path
from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for key in ('input','case','out'):parser.add_argument('--'+key,type=Path,required=True)
    parser.add_argument('--crest',type=Path,default=Path('/nfs/lambda_stor_01/homes/rzhu/miniforge3/envs/crest/bin/crest'))
    parser.add_argument('--threads',type=int,default=14);parser.add_argument('--max-hours',type=float,default=48)
    parser.add_argument('--energy-window-kcal',type=float,default=8.);parser.add_argument('--dry',action='store_true');args=parser.parse_args()
    if not 1<=args.threads<=64 or not 0<args.max_hours<=48 or args.energy_window_kcal<=0:parser.error('Invalid resource or sampling settings')
    with allocation_guard():
        capacity=host_capacity(args.threads,minimum_ram_gib=16);os.sched_setaffinity(0,cpu_allocation(args.threads))
    caps=dict(OMP_NUM_THREADS=str(args.threads),OMP_THREAD_LIMIT=str(args.threads),MKL_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',NUMEXPR_NUM_THREADS='1',XTB_THREADS='1')
    if any(os.environ.get(k)!=v for k,v in caps.items()):raise ValueError('Outer OpenMP and serial BLAS caps differ from the execution contract')
    from rdkit import Chem
    from redox.build.geometry import geometry_stereo_status
    from redox.validation.protonation import hydrogen_topology_status
    case=json.loads(args.case.read_text());mol=Chem.AddHs(Chem.MolFromSmiles(case['smiles']))
    charge=case['charge'];mult=case['mult'];electrons=sum(a.GetAtomicNum() for a in mol.GetAtoms())-charge
    if type(mult) is not int or mult<1 or (electrons-(mult-1))%2:raise ValueError('Invalid declared charge/spin parity')
    if hydrogen_topology_status(case['smiles'],args.input)!='ok' or geometry_stereo_status(case['smiles'],args.input) not in ('ok','not_applicable'):
        raise ValueError('Input graph/stereochemistry mismatch')
    args.out.mkdir(parents=True,exist_ok=False);(args.out/'input.xyz').write_bytes(args.input.read_bytes())
    xtb=args.crest.with_name('xtb');legacy=ROOT/'scripts/pipeline/sample_polymer_crest.py'
    spec=importlib.util.spec_from_file_location('frozen_crest_xyz_reader',legacy);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    sources=[Path(__file__).resolve(),legacy,args.input.resolve(),args.case.resolve(),args.crest.resolve(),xtb.resolve(),
        ROOT/'src/redox/workflows/resources.py',ROOT/'src/redox/build/geometry.py',ROOT/'src/redox/validation/protonation.py']
    hashes={str(p):sha(p) for p in sources};snapshot=args.out/'producer_snapshot';snapshot.mkdir()
    for i,p in enumerate(sources):
        if p.suffix in ('.py','.json','.xyz'):(snapshot/(str(i)+'_'+p.name)).write_bytes(p.read_bytes())
    argv=[str(args.crest),'input.xyz','--gfn2','--alpb','acetonitrile','--chrg',str(charge),'--uhf',str(mult-1),
        '--T',str(args.threads),'--ewin',str(args.energy_window_kcal),'--keepdir','--nosz','-xnam',str(xtb)]
    if args.dry:argv.append('--dry')
    status=dict(status='starting',host=os.uname().nodename,pid=os.getpid(),process_start=Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()[19],
        CPU_affinity=sorted(os.sched_getaffinity(0)),thread_caps=caps,capacity=capacity,case=case,argv=argv,source_sha256=hashes,
        started=time.time(),low_level_sampling_only=True,DFT_reranking_complete=False,ensemble_G_computed=False,
        MD_seed='Pinned CREST 3.0.2 default; independent realized trajectories retained. No documented iMTD seed override.')
    atomic_write_json(args.out/'status.json',status);child=None
    def stop(sig,frame):raise KeyboardInterrupt(str(sig))
    signal.signal(signal.SIGTERM,stop)
    try:
        with (args.out/'crest.log').open('w') as log:
            child=subprocess.Popen(argv,cwd=args.out,env=dict(os.environ,PATH=str(args.crest.parent)+os.pathsep+os.environ.get('PATH','')),stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
            status.update(status='running',predictor_pid=child.pid,predictor_start=Path(f'/proc/{child.pid}/stat').read_text().rsplit(')',1)[1].split()[19]);atomic_write_json(args.out/'status.json',status)
            code=child.wait(timeout=args.max_hours*3600)
        if code:raise RuntimeError('CREST terminated without successful completion')
        if any(sha(p)!=h for p,h in hashes.items()):raise ValueError('Sampling source/input changed')
        if args.dry:status['status']='dry_settings_check_complete_no_sampling'
        else:
            ensemble=args.out/'crest_conformers.xyz'
            if not ensemble.is_file():raise ValueError('No final conformer ensemble')
            checked=args.out/'checked_conformers';checked.mkdir();records=[]
            for i,frame in enumerate(module.xyz_frames(ensemble)):
                path=checked/f'conf_{i:05d}.xyz';path.write_text(frame)
                graph=hydrogen_topology_status(case['smiles'],path);stereo=geometry_stereo_status(case['smiles'],path)
                records.append(dict(index=i,path=str(path),sha256=sha(path),graph_QC=graph,stereo_QC=stereo,eligible_for_DFT_reranking=graph=='ok' and stereo in ('ok','not_applicable')))
            passed=sum(r['eligible_for_DFT_reranking'] for r in records)
            atomic_write_json(args.out/'conformers.json',dict(records=records,ensemble_sha256=sha(ensemble)))
            if not passed:raise ValueError('No conformers preserve the declared graph and stereochemistry')
            status.update(status='low_level_ensemble_complete_DFT_pending',conformers=len(records),topology_and_stereo_pass=passed)
    except BaseException as error:
        status.update(status='failed_or_interrupted_preserving_outputs',error=repr(error));raise
    finally:
        if child is not None and child.poll() is None:
            os.killpg(child.pid,signal.SIGTERM)
            try:child.wait(timeout=30)
            except subprocess.TimeoutExpired:os.killpg(child.pid,signal.SIGKILL);child.wait()
        status['finished']=time.time();atomic_write_json(args.out/'status.json',status)


if __name__=='__main__':main()
