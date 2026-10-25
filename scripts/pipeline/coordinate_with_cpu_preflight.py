#!/usr/bin/env python
"""Resume an existing controller with a CPU-safe dependency-probe adapter.

Physical producer files, plans, checkpoints and GPU validation remain unchanged.
The adapter is process-local; it does not edit the canonical fleet skill or
replace provenance files used by running calculations.
"""
import argparse
import inspect
from pathlib import Path
import sys
import time

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api


def main():
    parser=argparse.ArgumentParser(description=__doc__,add_help=False)
    parser.add_argument('--controller',required=True,choices=('coordinate_geometry_batch.py',
        'coordinate_smooth_geometry.py','coordinate_smooth_electronic_batch.py'))
    args,remaining=parser.parse_known_args()
    if '--batch' not in remaining:raise ValueError('An exact existing batch is required')
    batch=Path(remaining[remaining.index('--batch')+1]).resolve()
    upstream=ROOT/'scripts/pipeline'/args.controller
    helper=ROOT/'scripts/fleet/preflight_cpu_dependencies.py'
    proof={str(p):api.digest(p) for p in (Path(__file__).resolve(),upstream,helper)}
    original_loader,original_verify=api.load_module,api.verify_sources
    def verify(mapping):
        original_verify(proof);return original_verify(mapping)
    def load(name,path):
        verify({})
        module=original_loader(name,path)
        if Path(path).resolve()==ROOT/'scripts/pipeline/coordinate_geometry_batch.py':
            module.PREFLIGHT=helper
        elif Path(path).resolve()==ROOT/'scripts/pipeline/dispatch_smooth_energy_pilot.py':
            body=inspect.getsource(module.preflight)
            anchor="helper = ROOT.parent / '.codex-personal/skills/lambda-fleet/scripts/preflight_node.py'"
            if body.count(anchor)!=1:raise ValueError('Dependency-helper template changed')
            body=body.replace(anchor,'helper = Path('+repr(str(helper))+')',1)
            local={};exec(compile(body,str(Path(__file__).resolve()),'exec'),dict(vars(module)),local)
            module.preflight=local['preflight']
        return module
    api.verify_sources,api.load_module=verify,load
    history=batch/'CPU_preflight_transition';history.mkdir(exist_ok=True)
    atomic_write_json(history/(str(time.time_ns())+'.json'),dict(
        status='CPU_dependency_probe_adapter_selected',source_sha256=proof,
        upstream_controller=str(upstream),GPU_imports_deferred_to_reserved_workers=True,
        scientific_producer_bytes_changed=False,workers_or_GPU_reservations_signaled=False))
    module=load('source_bound_controller_CPU_preflight_transition',upstream)
    sys.argv=[str(upstream),*remaining]
    module.main()


if __name__=='__main__':main()
