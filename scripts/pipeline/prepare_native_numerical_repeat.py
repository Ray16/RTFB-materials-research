#!/usr/bin/env python3
"""Prepare one native repeat after a source-bound gradient reproducibility control."""
import argparse
import importlib.util
import json
import math
import os
from pathlib import Path
import re
import shutil

from redox.core.common import ROOT
from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json
from redox.workflows.archived_native_namespace import verify_spec
from redox.workflows.derived_native_namespace import prepare
from redox.workflows.submission_refusal import producer_inventory


def validate_control(analysis):
    values=[analysis['original_failed_relative_error'],analysis['repeated_relative_error_against_source_Hessian']]
    values.extend(v for r in analysis['points'] for v in (r['energy_difference_Ha'],r['maximum_gradient_difference_Ha_Bohr']))
    if any(type(v) not in (int,float) or not math.isfinite(v) for v in values):raise ValueError('Finite diagnostic evidence required')
    if (analysis['status']!='three_SCF_gradient_diagnostic_reconstructed_no_native_admission'
            or analysis['full_native_G_qualified'] or analysis['original_failed_relative_error']<=.01
            or analysis['repeated_relative_error_against_source_Hessian']>.01
            or len(analysis['points'])!=3 or max(abs(r['energy_difference_Ha']) for r in analysis['points'])>1e-9
            or max(r['maximum_gradient_difference_Ha_Bohr'] for r in analysis['points'])<=1e-8):
        raise ValueError('A completed same-energy, nonreproducing-gradient control is required')


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--recipe',type=Path,required=True);a=p.parse_args()
    from redox.workflows.resources import allocation_guard,host_capacity,cpu_allocation
    with allocation_guard():host_capacity(1,minimum_ram_gib=4);os.sched_setaffinity(0,cpu_allocation(1))
    config=json.loads(a.recipe.read_text());analysis_path=Path(config['diagnostic_analysis']);analysis=json.loads(analysis_path.read_text())
    validate_control(analysis)
    if any(sha256_file(k)!=v for k,v in analysis['source_sha256'].items()):raise ValueError('Diagnostic source changed')
    namespace=Path(analysis['original_namespace']);source=verify_spec(namespace);case=Path(source['attempt'])
    stage=json.loads((case/'experimental_native/status.json').read_text());old=json.loads((case/'plan.json').read_text())
    if stage['status']!='experimental_native_failed_preserved' or stage.get('error')!="ValueError('Alternate native sampled derivative gate failed')":
        raise ValueError('Only the diagnosed source derivative failure may enter this bounded repeat')
    definition=importlib.util.spec_from_file_location('numerical_repeat_owner_check',ROOT/'scripts/pipeline/observe_archived_native_namespace.py');owners=importlib.util.module_from_spec(definition);definition.loader.exec_module(owners)
    if owners.alive(stage['owner']) is not False:raise ValueError('Source native owner is live or unverified')
    policy=(ROOT.parent/'.config/gpu_reserve/config.env').read_text();hosts=re.search(r'^ALLOWED_HOSTS="([^"]+)"',policy,re.M)[1].split()
    inventory=producer_inventory(source,namespace,hosts)
    if set(inventory)!=set(hosts) or any(inventory.values()):raise ValueError('Matching native producer or unresolved fleet ownership')
    destination=Path(config['destination']);result=prepare(Path(source['source_energy_namespace']),destination)
    new=json.loads((Path(result['attempt'])/'plan.json').read_text())
    fields=['tasks','quantum_settings','experimental_model','experimental_method_id','experimental_native_settings',
        'scored_geometry_sha256','maximum_host_RSS_GiB','maximum_peak_GPU_MiB','execution_driver']
    if any(new[k]!=old[k] for k in fields):raise ValueError('Native repeat changed an exact scientific/state/geometry setting')
    proof=dict(status='one_full_native_repeat_prepared_after_gradient_nonreproducibility',
        original_namespace=str(namespace),original_namespace_sha256=sha256_file(namespace),
        diagnostic_analysis=str(analysis_path),diagnostic_analysis_sha256=sha256_file(analysis_path),
        producer_inventory=inventory,scientific_fields_compared=fields,scientific_fields_changed=False,
        full_native_G_qualified=False,automatic_retries=0,original_failed_outputs_preserved=True)
    atomic_write_json(destination/'numerical_repeat.json',proof);shutil.copyfile(__file__,destination/'repeat_preparer.py')
    print(json.dumps(dict(status=proof['status'],namespace=str(destination/'namespace.json'))))


if __name__=='__main__':main()
