"""Prepare a fresh ionic validation namespace from a proven closed full Hessian."""
import argparse
import json
import os
from pathlib import Path
import shutil

from redox.core.common import ROOT
from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json
from redox.validation.charged_hessian_source import source_evidence
from redox.workflows import molecular_native_namespace as execution


def prepare(namespace, output):
    owners=execution.load(ROOT/'scripts/pipeline/observe_archived_native_namespace.py','ionic_matrix_source_owners')
    spec,old,proof=source_evidence(namespace,owners.alive)
    output=Path(output).resolve()
    if output.exists() or any(output.is_relative_to(Path(spec[k])) for k in ('prior_root','archive')):
        raise FileExistsError('Preserve original evidence; use a new durable namespace')
    case=output/'case';sources=output/'sources';case.mkdir(parents=True);sources.mkdir()
    for name in ('dft','locks','repairs'):(case/name).mkdir()
    for path in Path(spec['sources']).iterdir():
        if path.is_file():shutil.copyfile(path,sources/path.name)
    replacements=[(ROOT/'scripts/validation/run_saved_polymer_hessian.py','worker.py'),
        (ROOT/'scripts/validation/run_screening_polymer_hessian.py','scientific_worker.py'),
        (ROOT/'src/redox/qm/saved_hessian.py','saved_matrix_helper.py'),
        (ROOT/'src/redox/validation/charged_hessian_source.py','source_evidence.py'),
        (Path(__file__),'recovery_preparer.py')]
    for path,name in replacements:shutil.copyfile(path,sources/name)
    atomic_write_json(sources/'source_matrix_evidence.json',proof)
    plan=dict(old);plan['source_sha256']=dict(old['source_sha256'])
    for key in ('resume_partial_Hessian','resume_components','resume_response'):
        plan.pop(key,None)
    origin=Path(spec['case'])
    saved=dict(namespace=str(namespace),source_case=str(origin),matrix=str(origin/'Hessian.npy'),
        center=str(origin/'points/0000/evaluation.json'),solver=str(origin/'execution/response_solver.json'),
        charge=old['charge'],mult=old['mult'],source_response_independently_checked=True,
        source_sha256=proof['source_sha256'])
    plan.update(Hessian_validation_continuation=True,charged_saved_Hessian_validation=True,
        saved_full_Hessian=saved,
        checkpoint_policy='Retain exact complete source Hessian; fresh center and every declared displaced-gradient, mode and partition check in a separate namespace',
        prior_failure_preserved=True,full_response_recomputed=False,native_G_qualified=False)
    plan['source_sha256'].update(proof['source_sha256'])
    plan['source_sha256'].update({str(p):sha256_file(p) for p in sources.iterdir() if p.is_file()})
    atomic_write_json(case/'plan.json',plan)
    new=dict(spec,case=str(case),sources=str(sources),plan_sha256=sha256_file(case/'plan.json'))
    atomic_write_json(output/'namespace.json',new);execution.verify(output/'namespace.json')
    return dict(status='closed_ionic_Hessian_fresh_validation_prepared',namespace=str(output/'namespace.json'),
        source_evidence=proof,quantum_jobs_submitted=0,native_G_qualified=False)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',type=Path,required=True);parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args()
    from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity
    with allocation_guard():host_capacity(1,minimum_ram_gib=4);os.sched_setaffinity(0,cpu_allocation(1))
    result=prepare(args.source.resolve(),args.out.resolve());print(json.dumps(result))
