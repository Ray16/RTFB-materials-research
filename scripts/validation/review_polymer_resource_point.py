#!/usr/bin/env python
"""Independently check completed storage/SCF pilots without promoting minima."""
import argparse
import json
import math
from pathlib import Path
import time

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.core.protocol import geom_sha1
from redox.workflows import polymer_campaign as api


def checked_numbers(data):
    result=data['result']
    values={key:result[key] for key in ('e_smd_eV','e_gas_eV')}
    values['gradient_max_Ha_Bohr']=data['gradient_max_Ha_Bohr']
    if any(type(v) not in (int,float) or not math.isfinite(v) for v in values.values()):
        raise ValueError('Each energy and gradient must be finite before aggregation')
    if any(result.get(key) is not True for key in ('converged_smd','converged_gas')) or values['gradient_max_Ha_Bohr']<0:
        raise ValueError('Both SCFs and a valid maximum gradient are required')
    if (result['charge'],result['mult'],result['solvent'])!=(0,1,'acetonitrile'):
        raise ValueError('Unexpected electronic state or medium')
    return values


def point_geometry(data):
    path=Path(data['result']['xyz'])
    expected=data['result'].get('energy_geom_sha1')
    if not expected or geom_sha1(path)!=expected:
        raise ValueError('The energy input geometry differs or is missing')
    return path


def review(campaign,out):
    plan=api.read_plan(campaign);api.verify_sources(plan['producer_sha256']);api.verify_sources(plan['input_sha256'])
    point=campaign/'energy/pilot/energy_and_gradient.json';status_path=campaign/'energy/pilot/status.json'
    status=json.loads(status_path.read_text());data=json.loads(point.read_text());values=checked_numbers(data)
    if status['status']!='memory_pilot_passed' or not 0<status['observed_peak_GPU_MiB']<=plan['maximum_memory_pilot_MiB']:
        raise ValueError('The normal pilot/resource gates must pass')
    api.verify_sources(status['output_sha256'])
    geometry=point_geometry(data)
    input_path=campaign/'live_inputs'/plan['tasks'][0]['key']/'input.xyz'
    if api.digest(geometry)!=api.digest(input_path):raise ValueError('The point input belongs to another task')
    guard_path=campaign/'host_guard/status.json';guard=json.loads(guard_path.read_text())
    if (guard['status']!='program_returned_review_required' or guard.get('error')
            or not 0<guard['peak_RSS_bytes']<=guard['maximum_RSS_GiB']*2**30
            or guard['minimum_free_GiB']!=32
            or not guard['started']<=status['started']<=status['finished']<=guard['finished']):
        raise ValueError('Host-resource completion is missing or invalid')
    owner=json.loads((campaign/'residency_guard/owner.json').read_text())
    if owner['pid']!=status['pid'] or api.load_module('resource_review_owner',ROOT/'scripts/pipeline/run_polymer_spacing.py').owner_alive(owner):
        raise ValueError('The exact scientific owner must have ended')
    sources={str(p):api.digest(p) for p in (campaign/'plan.json',point,status_path,guard_path,Path(__file__).resolve())}
    differences=None
    if plan.get('comparison_reference_profile'):
        reference=Path(plan['comparison_reference_profile'])/'energy/pilot/energy_and_gradient.json'
        previous=json.loads(reference.read_text());reference_values=checked_numbers(previous)
        if api.digest(geometry)!=api.digest(point_geometry(previous)):raise ValueError('Reference coordinates differ')
        differences={key:abs(value-reference_values[key]) for key,value in values.items()}
        if any(v>(1e-6 if key=='gradient_max_Ha_Bohr' else 1e-5) for key,v in differences.items()):
            raise ValueError('The fixed-point energies/maximum force differ from the preserved reference')
        sources[str(reference)]=api.digest(reference)
    out.mkdir(parents=True,exist_ok=False)
    report=dict(status='electronic_resource_point_independently_reviewed',checked=time.time(),
        source=str(campaign),task_key=plan['tasks'][0]['key'],values=values,
        larger_device_reference_differences=differences,observed_peak_GPU_MiB=status['observed_peak_GPU_MiB'],
        peak_host_RSS_bytes=guard['peak_RSS_bytes'],host_RSS_bound_GiB=guard['maximum_RSS_GiB'],
        source_sha256=sources,optimized_geometry_qualified=False,native_Hessian_qualified=False,G_computed=False,
        scope='Fixed-input SCF/resource evidence; reference gradient comparison is a scalar maximum, not a full vector')
    atomic_write_json(out/'summary.json',report)
    print(json.dumps({key:report[key] for key in ('status','observed_peak_GPU_MiB','larger_device_reference_differences')}))


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--campaign',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();review(a.campaign.resolve(),a.out.resolve())


if __name__=='__main__':main()
