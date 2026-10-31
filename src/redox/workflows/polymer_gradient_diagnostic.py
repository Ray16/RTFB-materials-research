"""Source-bound three-SCF gradient reproduction, independent of native admission."""
import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess

import numpy as np

from redox.core.common import ROOT
from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json
from redox.validation.namespace_sources import NamespaceSources
from redox.workflows.archived_native_namespace import verify_spec
from redox.workflows import molecular_native_namespace as execution


def select_points(raw):
    checks=json.loads((raw/'derivative_checks.json').read_text())['checks']
    failed=max(checks,key=lambda r:r['relative_error'])
    if failed['relative_error']<=.01:raise ValueError('The requested diagnostic must have a failed saved column')
    paths=sorted((raw/'SCF_evaluations').glob('*/evaluation.json'));base=json.loads(paths[0].read_text());xyz=np.asarray(base['coordinates_Bohr'])
    selected=[paths[0]]
    for sign in (-1,1):
        target=xyz.copy();target[failed['atom'],failed['axis']]+=sign*failed['step_Bohr']
        matches=[p for p in paths if np.max(abs(np.asarray(json.loads(p.read_text())['coordinates_Bohr'])-target))<1e-9]
        if len(matches)!=1:raise ValueError('One exact saved displaced SCF required')
        selected.append(matches[0])
    return failed,selected


def prepare(namespace,output):
    if output.exists() or output.resolve().is_relative_to(ROOT):raise FileExistsError('Use a fresh durable diagnostic destination')
    source=verify_spec(namespace);reader=NamespaceSources(source);original=Path(source['attempt']);parent=json.loads((original/'plan.json').read_text())
    if any(output.resolve().is_relative_to(Path(source[key]).resolve()) for key in ('archive','attempt','raw_parent')):
        raise ValueError('Diagnostic outputs must remain outside every scientific source tree')
    stage=json.loads((original/'experimental_native/status.json').read_text())
    if stage['status']!='experimental_native_failed_preserved' or 'finished' not in stage:raise ValueError('A closed failed native record is required')
    if len(parent['tasks'])!=1 or (parent['tasks'][0]['charge'],parent['tasks'][0]['mult'])!=(0,1):raise ValueError('Declared neutral singlet required')
    reader.verify_all(parent['producer_sha256']);reader.verify_all(parent['input_sha256'])
    raw=reader.path(source['original_raw_root']);failed,points=select_points(raw)
    output.mkdir(parents=True);case=output/'case';sources=output/'sources';case.mkdir();sources.mkdir()
    for name in ('dft','locks','repairs'):(case/name).mkdir()
    for path,name in [(ROOT/'scripts/validation/run_polymer_gradient_diagnostic.py','worker.py'),
        (ROOT/'src/redox/qm/scf_checkpoints.py','checkpoints.py'),
        (ROOT/'scripts/fleet/bounded_peak_host_python_worker.py','host_guard.py'),(Path(__file__),'preparer.py')]:shutil.copyfile(path,sources/name)
    inputs=[namespace,original/'plan.json',original/'experimental_native/status.json',raw/'geometry.xyz',raw/'derivative_checks.json',raw/'hessian.npy',*points,*sources.glob('*')]
    for name in ('scripts/fleet/exclusive_python_worker.py','scripts/fleet/tensor_runtime_worker.py','environments/cutensor22/runtime.json'):
        inputs.append(Path(source['archive'])/name)
    bindings={**{str(reader.path(p)):h for p,h in {**parent['producer_sha256'],**parent['input_sha256']}.items()},
        **{str(p):sha256_file(p) for p in inputs}}
    plan=dict(role='polymer_gradient_diagnostic',source_id='polymer_gradient_diagnostic',state='rest',
        original_namespace=str(namespace),original_namespace_sha256=sha256_file(namespace),failed_column=failed,
        symbols=[s.split()[0] for s in (raw/'geometry.xyz').read_text().splitlines()[2:]],
        points=[dict(evaluation=str(p),source_point=p.parent.name) for p in points],
        source_sha256=bindings,maximum_host_RSS_GiB=80,minimum_free_GiB=32,maximum_peak_GPU_MiB=31744,
        native_G_qualified=False,material_admission=False,original_scientific_records_writable=False)
    atomic_write_json(case/'plan.json',plan)
    spec=dict(status='polymer_gradient_diagnostic_prepared',prior_root=source['prior_root'],archive=source['archive'],
        case=str(case),sources=str(sources),source_id=plan['source_id'],state='rest',maximum_host_RSS_GiB=80,minimum_free_GiB=32,plan_sha256=sha256_file(case/'plan.json'))
    spec['scientific_sources_read_only']=True
    atomic_write_json(output/'namespace.json',spec);reader.recheck();return spec


def probe(path):
    spec=execution.verify(path);runtime=json.loads((Path(spec['archive'])/'environments/cutensor22/runtime.json').read_text())
    code=('import json,ast,pathlib,importlib.util,importlib.metadata as m;import ase,rdkit,numpy,psutil;'
        'wanted='+repr(runtime['versions'])+';assert {k:m.version(k) for k in wanted}==wanted;'
        'tree=ast.parse(pathlib.Path('+repr(spec['sources']+'/worker.py')+').read_text());'
        'modules=sorted({n.module for n in ast.walk(tree) if isinstance(n,ast.ImportFrom) and n.module and n.module.startswith("redox.")});'
        'assert all(importlib.util.find_spec(name) is not None for name in modules),modules;'
        'definition=importlib.util.spec_from_file_location("component_witness",'+repr(spec['prior_root']+'/scripts/pipeline/smooth_cds_control.py')+');'
        'model=importlib.util.module_from_spec(definition);definition.loader.exec_module(model);'
        'witness=model.component_evidence(pathlib.Path('+repr(spec['prior_root']+'/results/validation/smooth_CDS_2026-10-25/pilot')+'));'
        'print(json.dumps(dict(status="diagnostic_CPU_dependencies_checked",versions=wanted,worker_modules=modules,component_source_bindings=len(witness))))')
    run=subprocess.run(execution.command(spec)+[runtime['python'],'-c',code],capture_output=True,text=True,check=True,timeout=180)
    report=dict(status='polymer_gradient_diagnostic_node_sources_checked',host=os.uname().nodename.split('.')[0],
        namespace_sha256=sha256_file(path),checks=json.loads(run.stdout),GPU_started=False)
    atomic_write_json(path.with_name('node_probe.json'),report);return report


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('action',choices=('prepare','probe','launch'))
    p.add_argument('--spec',type=Path,required=True);p.add_argument('--out',type=Path);p.add_argument('--gpu',type=int);a=p.parse_args()
    if a.action!='launch':
        from redox.workflows.resources import allocation_guard,host_capacity,cpu_allocation
        with allocation_guard():host_capacity(1,minimum_ram_gib=4);os.sched_setaffinity(0,cpu_allocation(1))
    if a.action=='prepare':result=prepare(a.spec.resolve(),a.out.resolve())
    elif a.action=='probe':result=probe(a.spec.resolve())
    else:
        child,result=execution.launch(a.spec.resolve(),a.gpu,a.out.resolve());print(json.dumps(result),flush=True)
        result.update(returncode=child.wait(),status='polymer_gradient_diagnostic_returned_analysis_pending');atomic_write_json(a.out/'completion.json',result)
    print(json.dumps(result))
