#!/usr/bin/env python
"""Run existing thermochemistry with an explicit hash-qualified xTB prefix."""
import argparse
import json
import os
from pathlib import Path
import sys

from redox.qm.xtb_runtime import qualified_xtb_prefix


def probe(runtime):
    from redox.qm import qrrho
    masses=qrrho.atomic_masses([1,6,7,8,17]).tolist()
    thermal=qrrho.kernel_thermal([50.,100.,500.],[1.,.5,.3],2.0e-46,30.,1,cutoff_cm=50)
    worker=qrrho.kernel_execution();maps=Path('/proc')/str(worker['pid'])/'maps'
    if str(Path(runtime['library']).resolve()) not in maps.read_text():
        raise RuntimeError('Thermodynamic helper loaded a different library path')
    return dict(status='qualified_xtb_relocated_kernel_probe_checked',runtime=runtime,masses_amu=masses,
        thermal=thermal,worker=worker,loaded_library_path_verified=True,quantum_calculations_run=False,
        native_state_qualified=False,full_platform_portability_inferred=False)


def configure_thread_pools(mode,forward):
    parser=argparse.ArgumentParser(add_help=False)
    parser.add_argument('--threads',type=int,default=2)
    parsed,_=parser.parse_known_args(forward)
    threads=1 if mode=='probe' else parsed.threads
    if threads<1:raise ValueError('Positive thread count required')
    for key in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS','XTB_THREADS'):
        os.environ[key]=str(threads)
    return threads


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--xtb-prefix',type=Path,
        default=Path(os.environ['REDOX_XTB_PREFIX']) if os.environ.get('REDOX_XTB_PREFIX') else None)
    p.add_argument('--mode',choices=['probe','projected','native'],required=True);p.add_argument('--out',type=Path)
    args,forward=p.parse_known_args()
    if forward and forward[0]=='--':forward=forward[1:]
    if args.xtb_prefix is None:p.error('Provide --xtb-prefix or REDOX_XTB_PREFIX')
    threads=configure_thread_pools(args.mode,forward)
    if args.mode in ('probe','projected'):
        from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity
        with allocation_guard():host_capacity(threads);os.sched_setaffinity(0,cpu_allocation(threads))
    else:
        from redox.workflows.resources import verify_gate_owner
        verify_gate_owner()
        if len(os.sched_getaffinity(0))!=threads:p.error('Native mode requires the reserved, pinned CPU allocation')
    with qualified_xtb_prefix(args.xtb_prefix) as runtime:
        if args.mode=='probe':
            if forward:p.error('Probe mode does not accept additional arguments')
            result=probe(runtime)
            if args.out:
                if args.out.exists():raise FileExistsError('Preserve the existing probe')
                args.out.parent.mkdir(parents=True,exist_ok=True);args.out.write_text(json.dumps(result,indent=2)+'\n')
            print(json.dumps(result,indent=2))
        elif args.mode=='projected':
            from redox.qm import projected_thermal
            if args.out is not None:forward.extend(['--out',str(args.out)])
            projected_thermal.main(forward)
        else:
            from redox.qm import native_thermal
            if args.out is not None:p.error('--out is only used by probe or projected mode')
            sys.argv=[sys.argv[0],*forward];native_thermal.main()


if __name__=='__main__':main()
