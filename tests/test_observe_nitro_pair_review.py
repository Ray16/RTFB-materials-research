import importlib.util
import json
from pathlib import Path

import pytest


source = Path(__file__).resolve().parents[1] / 'scripts/validation/observe_nitro_pair_review.py'
spec = importlib.util.spec_from_file_location('nitro_pair_observer_test', source)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def data(tmp_path):
    identity = {'model':'m','checkpoint_sha256':'c','reference_sha256':{},'task_name':'omol'}
    case = {'id':'a','charge':0,'mult':1,'atoms':1,'frames':['one']}
    plan = dict(model_identity=identity, cases=[case], minimum_energy_tolerance_eV=.002,
        energy_tolerance_eV_per_atom=.0001, force_RMSE_tolerance_eV_A=.002,
        force_max_tolerance_eV_A=.01, input_sha256={}, producer_sha256={})
    plan_path = tmp_path/'plan.json'; plan_path.write_text(json.dumps(plan))
    frame = {'energy_eV':-2.,'forces_eV_A':[[0.,0.,0.]]}
    row = dict(id='a',charge=0,mult=1,atoms=1,frames=[frame],actual_mode='reference')
    common = dict(plan_sha256=module.api.digest(plan_path),host='node',GPU='1',capability=[9,0],
                  status='UMA_point_evaluations_complete_cross_mode_review_pending',model_identity=identity)
    ref = dict(common, mode='reference', explicit_reference_mode_checks=1,cases=[row])
    nitro_row = dict(row,actual_mode='nitro',execution={'actual_mode':'nitro','executor':'fairchem.uma_nitro.Model'})
    nitro = dict(common, mode='nitro',cases=[nitro_row])
    pair = tmp_path/'pair'; (pair/'reference').mkdir(parents=True); (pair/'nitro').mkdir()
    (pair/'reference/status.json').write_text(json.dumps(ref)); (pair/'nitro/status.json').write_text(json.dumps(nitro))
    guard = tmp_path/'guard'; guard.mkdir()
    (guard/'owner.json').write_text(json.dumps({'guard_stopped':True,'host':'node','GPU':1}))
    (guard/'state.json').write_text(json.dumps({'foreign_resident_PIDs':[],'own_residency_observed':True}))
    return plan_path,pair,guard


def test_fixed_pair_is_recomputed_without_claiming_optimizer_qualification(tmp_path):
    plan,pair,guard = data(tmp_path)
    result = module.independently_review(plan,pair,guard)
    assert result['status'] == 'independent_fixed_device_pair_passed'
    assert result['frames'] == 1 and not result['production_optimizer_qualified']


@pytest.mark.parametrize('mutation', ['wrong_gpu','not_nitro','foreign_gpu','unverified_reference',
                                       'wrong_guard','observer_error','diagnostic'])
def test_incomplete_or_mixed_evidence_is_rejected(tmp_path,mutation):
    plan,pair,guard = data(tmp_path)
    path = pair/'nitro/status.json'; record=json.loads(path.read_text())
    if mutation == 'wrong_gpu': record['GPU']='2'
    if mutation == 'not_nitro': record['cases'][0]['execution']['actual_mode']='other'
    if mutation in ('wrong_gpu','not_nitro'): path.write_text(json.dumps(record))
    if mutation == 'foreign_gpu': (guard/'state.json').write_text(json.dumps({'foreign_resident_PIDs':[3],'own_residency_observed':True}))
    if mutation == 'unverified_reference':
        path=pair/'reference/status.json';record=json.loads(path.read_text());record['explicit_reference_mode_checks']=0;path.write_text(json.dumps(record))
    if mutation == 'wrong_guard': (guard/'owner.json').write_text(json.dumps({'guard_stopped':True,'host':'node','GPU':9}))
    if mutation == 'observer_error': (guard/'observation_error.json').write_text('{}')
    if mutation == 'diagnostic':
        path=pair/'nitro/status.json';record=json.loads(path.read_text());record['diagnostic_only']=True;path.write_text(json.dumps(record))
    with pytest.raises(ValueError):
        module.independently_review(plan,pair,guard)
