#!/usr/bin/env python
"""Place one energy pilot with its actual RAM bound in CPU-only preflight.

Insufficient host capacity is contention before any process is spawned. The
physical worker, numerical settings and independently checked geometry stay fixed.
"""
import inspect
from pathlib import Path

from redox.core.common import ROOT
from redox.workflows import polymer_campaign as api


def capacity_placement(profile):
    base=api.load_module('energy_capacity_reserved_placement',ROOT/'scripts/pipeline/coordinate_geometry_batch.py')
    body=inspect.getsource(base.launch_local)
    anchor="        with (attempt / 'worker.log').open('x') as log:"
    capacity="        capacity = host_capacity(2, minimum_ram_gib=bound + 32)"
    if body.count(anchor)!=1 or body.count(capacity)!=1:raise ValueError('Capacity placement template changed')
    body=body.replace(capacity,"""        try:
            capacity = host_capacity(2, minimum_ram_gib=bound + 32)
        except RuntimeError as error:
            return dict(status='resource_capacity_unavailable', host=host, reason=str(error))""",1)
    inserted="""        if plan.get('runtime_profile') != str(_PROFILE):
            raise ValueError('Unexpected energy runtime profile')
        argv.insert(argv.index('env')+1, 'LD_LIBRARY_PATH='+_profile['library_path'])
        index=argv.index(manifest['worker'])
        argv[index:index]=[str(ROOT/'scripts/fleet/tensor_runtime_worker.py'), '--profile', str(_PROFILE),
            '--out', str(campaign/'tensor_runtime')]
"""
    namespace=dict(vars(base),_PROFILE=ROOT/'environments/cutensor22/runtime.json',_profile=profile);local={}
    exec(compile(body.replace(anchor,inserted+anchor,1),str(Path(__file__).resolve()),'exec'),namespace,local)
    base.launch_local=local['launch_local'];return base


def capacity_preflight(helper):
    body=inspect.getsource(helper.preflight)
    old="helper = ROOT.parent / '.codex-personal/skills/lambda-fleet/scripts/preflight_node.py'"
    if body.count(old)!=1 or body.count("'65536'")!=1:raise ValueError('RAM preflight template changed')
    body=body.replace(old,"helper = ROOT / 'scripts/fleet/preflight_cpu_dependencies.py'",1)
    body=body.replace("    path = batch / (host + '_preflight.json')","    required = str(int((json.loads((batch/'batch.json').read_text())['maximum_RSS_GiB']+32)*1024))\n    path = batch / (host + '_preflight.json')",1)
    body=body.replace("if path.exists() and time.time() - path.stat().st_mtime < 600:","if path.exists() and time.time() - path.stat().st_mtime < 600 and json.loads(path.read_text()).get('ok'):",1)
    body=body.replace("return data['ok'] and data['python'] == profile['python']","return data['ok'] and data['python'] == profile['python'] and data.get('minimum_RAM_MiB') == int(required)",1)
    body=body.replace("'65536'",'required',1)
    body=body.replace("data = json.loads(result.stdout); atomic_write_json(path, data)","data = json.loads(result.stdout); data['minimum_RAM_MiB']=int(required); atomic_write_json(path, data)",1)
    namespace={};exec(compile(body,str(Path(__file__).resolve()),'exec'),dict(vars(helper)),namespace)
    return namespace['preflight']


def main():
    helper=api.load_module('capacity_energy_pilot_contract',ROOT/'scripts/pipeline/dispatch_smooth_energy_pilot.py')
    helper.SCRIPT=Path(__file__).resolve();helper.WORKER=ROOT/'scripts/pipeline/smooth_polymer_energy_bounded.py'
    helper.placement_module=capacity_placement;helper.preflight=capacity_preflight(helper);helper.main()


if __name__=='__main__':main()
