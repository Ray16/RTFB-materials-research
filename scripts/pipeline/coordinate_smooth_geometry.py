#!/usr/bin/env python
"""Use the qualified queue/resource contract with an explicit model/runtime."""
import inspect
import json
from pathlib import Path
import sys

from redox.core.common import ROOT
from redox.workflows import polymer_campaign as api

PROFILE=ROOT/'environments/cutensor22/runtime.json'


def main():
    profile=json.loads(PROFILE.read_text())
    if Path(sys.executable)!=Path(profile['python']):
        raise ValueError('Use the declared isolated runtime interpreter for this coordinator')
    batch=Path(sys.argv[sys.argv.index('--batch')+1]).resolve()
    manifest=json.loads((batch/'batch.json').read_text())
    if (manifest['stage']!='geometry' or manifest['worker']!=str(ROOT/'scripts/pipeline/run_smooth_geometry.py')
            or manifest['primary_candidate_admission']):
        raise ValueError('Only the declared alternate-model geometry batch is supported')
    source=ROOT/'scripts/pipeline/coordinate_geometry_batch.py'
    module=api.load_module('unchanged_model_geometry_queue_contract',source)
    module.SCRIPT=Path(__file__).resolve()
    # Qualify imports using the actual loader path on each remote node.
    body=inspect.getsource(module.preflight)
    anchor="    result = subprocess.run(argv, capture_output=True, text=True, timeout=60)"
    if body.count(anchor)!=1:raise ValueError('Node preflight template changed')
    body=body.replace(anchor,"    argv += ['--library-path', _profile['library_path']]\n"+anchor,1)
    body=body.replace("'49152'","'65536'")
    namespace=dict(vars(module),_profile=profile);local={};exec(compile(body,str(Path(__file__).resolve()),'exec'),namespace,local)
    module.preflight=local['preflight']
    # The same-process guards must wrap the actual contraction-runtime worker.
    body=inspect.getsource(module.launch_local)
    anchor="        with (attempt / 'worker.log').open('x') as log:"
    if body.count(anchor)!=1:raise ValueError('Atomic geometry placement template changed')
    inserted="""        if plan.get('runtime_profile') != str(_PROFILE):
            raise ValueError('Task runtime differs from the model batch')
        argv.insert(argv.index('env')+1, 'LD_LIBRARY_PATH='+_profile['library_path'])
        index=argv.index(manifest['worker'])
        argv[index:index]=[str(ROOT/'scripts/fleet/tensor_runtime_worker.py'), '--profile', str(_PROFILE),
            '--out', str(campaign/'tensor_runtime')]
"""
    body=body.replace(anchor,inserted+anchor,1)
    namespace=dict(vars(module),_profile=profile,_PROFILE=PROFILE);local={}
    exec(compile(body,str(Path(__file__).resolve()),'exec'),namespace,local);module.launch_local=local['launch_local']
    body=inspect.getsource(module.coordinate)
    old='stationary_geometry_checked_active_parent_and_native_minimum_pending'
    if body.count(old)!=1:raise ValueError('Geometry-result status contract changed')
    body=body.replace(old,'experimental_smooth_geometry_checked_energy_and_native_QC_pending',1)
    namespace=dict(vars(module));local={};exec(compile(body,str(Path(__file__).resolve()),'exec'),namespace,local)
    module.coordinate=local['coordinate']
    module.main()


if __name__=='__main__':main()
