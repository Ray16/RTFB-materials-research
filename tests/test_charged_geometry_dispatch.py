"""Storage admission and ownership failures must preserve scientific attempts."""
import importlib.util
import json
from pathlib import Path

import pytest


@pytest.fixture
def module():
    path=Path(__file__).resolve().parents[1]/'scripts/fleet/run_charged_polymer_geometry_panel.py'
    spec=importlib.util.spec_from_file_location('test_ionic_geometry_dispatch',path)
    result=importlib.util.module_from_spec(spec);spec.loader.exec_module(result);return result


def case(tmp_path,name,status=None):
    root=tmp_path/name;folder=root/'case';folder.mkdir(parents=True)
    namespace=root/'namespace.json';namespace.write_text(json.dumps(dict(case=str(folder))))
    (folder/'plan.json').write_text(json.dumps(dict(maximum_durable_output_bytes=2**30)))
    if status is not None:(folder/'status.json').write_text(json.dumps(status))
    return dict(namespace=str(namespace),key=name)


def test_aggregate_storage_reserves_unknown_owners_and_deduplicates(module,tmp_path,monkeypatch):
    row=case(tmp_path,'a',dict(status='running',owner=dict(pid=1)))
    monkeypatch.setattr(module,'verify',lambda p:json.loads(p.read_text()))
    record=module.storage_headroom([row],[row],2**30,lambda owner:None,34*2**30)
    assert 0<record['remaining_existing_budget_bytes']<=2**30
    assert record['sufficient']
    assert not module.storage_headroom([row],[row],2**30,lambda owner:None,33*2**30)['sufficient']


def test_closed_failed_attempt_does_not_reserve_future_writes(module,tmp_path,monkeypatch):
    row=case(tmp_path,'closed',dict(status='failed',finished=1,owner=dict(pid=1)))
    monkeypatch.setattr(module,'verify',lambda p:json.loads(p.read_text()))
    assert module.storage_headroom([row],[],0,lambda owner:False,32*2**30)['sufficient']


def test_disappeared_scientific_owner_is_preserved(module,tmp_path):
    row=case(tmp_path,'gone',dict(status='running',owner=dict(pid=1)))
    assert module.inspect(row,lambda owner:False)[0]=='stopped_preserved'
    assert module.inspect(row,lambda owner:None)[0]=='ownership_unverified'


def test_uncertain_or_dead_submission_is_never_relaunched(module,tmp_path):
    row=case(tmp_path,'uncertain');claim=Path(row['namespace']).parent/'dispatch_claim';claim.mkdir()
    assert module.inspect(row,lambda owner:False)[0]=='uncertain_submission_preserved'
    (claim/'try_001.json').write_text(json.dumps(dict(owner=dict(pid=1))))
    assert module.inspect(row,lambda owner:False)[0]=='closed_submission_preserved'


@pytest.mark.parametrize('change',['PID_reused','foreign_GPU','overlapping_CPU','missing_thread_cap','wrong_worker'])
def test_startup_verification_requires_actual_process_evidence(change):
    path=Path(__file__).resolve().parents[1]/'scripts/validation/observe_charged_geometry_panel.py'
    spec=importlib.util.spec_from_file_location('ionic_start_checks',path)
    observer=importlib.util.module_from_spec(spec);spec.loader.exec_module(observer)
    row=dict(stage=dict(owner=dict(pid=7,process_start='123'),status='charged_polymer_solution_optimization_running'),
        launch=dict(device=dict(uuid='gpu-a'),CPU_affinity=[2,3]),worker='/worker.py',plan='/plan.json',SCF_progress=dict(cycle=1))
    process=dict(process_start='123',cpus=[2,3],argv=['python','/worker.py','--plan','/plan.json'],
        env={k:'2' for k in ['OMP_NUM_THREADS','OMP_THREAD_LIMIT','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS','XTB_THREADS']})
    process['env']['NUMPY_MADVISE_HUGEPAGE']='0';resident=[dict(uuid='gpu-a',pid=7)]
    assert all(observer.checks(row,process,resident).values())
    if change=='PID_reused':process['process_start']='456'
    if change=='foreign_GPU':resident.append(dict(uuid='gpu-a',pid=9))
    if change=='overlapping_CPU':process['cpus']=[2,4]
    if change=='missing_thread_cap':del process['env']['OMP_THREAD_LIMIT']
    if change=='wrong_worker':process['argv']=['python','/unrelated.py']
    assert not all(observer.checks(row,process,resident).values())
