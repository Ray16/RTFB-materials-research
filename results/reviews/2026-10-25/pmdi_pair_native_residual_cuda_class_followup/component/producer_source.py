#!/usr/bin/env python
"""Reconstruct completed polymer G components without inferring a complete cycle."""
import argparse
import json
import math
import os
from pathlib import Path
import time

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.core.protocol import load_record,record_path
from redox.qm.native_thermal import cache_path,load_native_thermal
from redox.workflows import polymer_campaign as api
from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity


def review(campaign,out):
    with allocation_guard():host_capacity(1);os.sched_setaffinity(0,cpu_allocation(1))
    plan=api.read_plan(campaign);api.verify_sources(plan['producer_sha256']);api.verify_sources(plan['input_sha256'])
    records,invalid=api.admitted_records(campaign,'thermal')
    if invalid or not records:raise ValueError('Completed normally admitted native components required')
    equations_path=ROOT/'scripts/validation/audit/native_partition.py'
    equations=api.load_module('independent_component_nuclear_equations',equations_path)
    geometry_path=ROOT/'scripts/validation/audit_polymer_geometries.py'
    geometry=api.load_module('independent_component_geometry',geometry_path)
    sources={str(p):api.digest(p) for p in (Path(__file__).resolve(),campaign/'plan.json',equations_path,geometry_path)}
    rows=[]
    for key,record in records.items():
        task=record['task'];state=Path(record['state_dir']);q,m=task['charge'],task['mult']
        energy=load_record(state,q,m);native=load_native_thermal(state,q,m)
        surface,_=cache_path(state,q,m);partition=equations.read_partition(surface,state/'opt.xyz')
        if energy is None or native is None or partition is None:
            raise ValueError('Missing admitted or independently reconstructed parent')
        audit=geometry.audit(task,state/'opt.xyz')
        if audit['issues']:raise ValueError('Component graph/geometry audit failed: '+repr(audit['issues']))
        G=energy['e_smd_eV']+partition['g_thermal_total_eV']
        if not math.isfinite(G) or abs(G-record['G_eV'])>1e-8:
            raise ValueError('Native component free energy does not reconstruct')
        receipt=campaign/'thermal/done'/(key+'.json')
        status=Path(json.loads(receipt.read_text())['status_path'])
        for p in (receipt,status,surface,record_path(state,q,m),state/'opt.xyz'):sources[str(p)]=api.digest(p)
        sources.update(record['output_sha256']);sources.update(partition['evidence_sha256'])
        rows.append(dict(key=key,occupied_sites=task.get('occupied_sites'),state_dir=str(state),
            E_eV=energy['e_smd_eV'],G_eV=G,
            G_by_rotor_cutoff_eV={str(k):energy['e_smd_eV']+partition['nuclear_by_cutoff_eV'][str(k)]+native['g_spin_partition_eV'] for k in (25,50,100)},
            internal_modes=len(native['frequencies_cm']),minimum_frequency_cm=min(native['frequencies_cm']),
            gradient_max_Ha_Bohr=native['gradient_max'],maximum_Hessian_asymmetry=native['maximum_asymmetry'],
            maximum_sampled_column_relative_error=max(c['relative_error'] for c in native['derivative_checks']),geometry_audit=audit))
    api.verify_sources(sources)
    out.mkdir(parents=True,exist_ok=False)
    result=dict(status='available_exact_polymer_native_components_independently_reconstructed',checked=time.time(),
        checked_components=len(rows),expected_campaign_components=len(plan['tasks']),states=rows,source_sha256=sources,
        matched_cycle_established=False,conformer_convergence_established=False,QM_recomputed=False,
        normal_admission_rechecked=True,scope='Exact state energy, independent Hessian projection/nuclear equations and geometry; not a complete occupancy or selection decision')
    atomic_write_json(out/'summary.json',result);(out/'producer_source.py').write_bytes(Path(__file__).read_bytes())
    print(json.dumps(dict(status=result['status'],checked=len(rows),G_eV=[r['G_eV'] for r in rows])))


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--campaign',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();review(a.campaign.resolve(),a.out.resolve())


if __name__=='__main__':main()
