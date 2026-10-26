import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace

import pytest


ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('review_uma_pilot_test', ROOT/'scripts/validation/review_external_uma_pilot.py')
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def fixture(tmp_path, monkeypatch):
    campaign = tmp_path/'campaign'; attempt = tmp_path/'attempt'
    for p in [campaign/'output/nitro/g/neu',campaign/'completed',campaign/'lib/g/conformers',
              attempt/'residency_guard',attempt/'host_guard',attempt/'gpu_memory']:
        p.mkdir(parents=True)
    (campaign/'plan.json').write_text('{}')
    policy_path = ROOT/'config/uma_execution.json'; policy=json.loads(policy_path.read_text())
    result = dict(model_identity={'checkpoint_sha256':policy['checkpoint_sha256'],'reference_sha256':policy['reference_sha256']},
        UMA_execution={'actual_mode':'nitro','policy_sha256':module.api.digest(policy_path)},energy_eV=-1.,fmax=.01)
    rp=campaign/'output/nitro/g/neu/result.json';rp.write_text(json.dumps(result))
    (campaign/'output/nitro/g/neu/relaxed.xyz').write_text('geometry')
    (campaign/'completed/neu.json').write_text(json.dumps({'status':'UMA_control_worker_returned_review_required','result_sha256':module.api.digest(rp),'owner':{'pid':123,'process_start':'start'}}))
    (attempt/'residency_guard/owner.json').write_text(json.dumps({'guard_stopped':True,'pid':123,'GPU':'1','process_start':'start','CPU_affinity':[1,2]}))
    (attempt/'residency_guard/state.json').write_text(json.dumps({'foreign_resident_PIDs':[],'own_residency_observed':True}))
    (attempt/'host_guard/status.json').write_text(json.dumps({'status':'program_returned_review_required','peak_RSS_bytes':1024,'pid':123,'CPU_affinity':[1,2]}))
    (attempt/'gpu_memory/status.json').write_text(json.dumps({'status':'GPU_memory_observation_program_returned_review_required',
        'own_samples':2,'sampling_errors':[],'peak_process_memory_MiB':50,'pid':123,'GPU':'1'}))
    plan = {'group':'g','library':str(campaign/'lib')}
    monkeypatch.setattr(module.api,'load_module',lambda *a:SimpleNamespace(validate_campaign=lambda c:(plan,[{'state':'neu'}])))
    monkeypatch.setattr(module.uma,'cache_errors',lambda *a:[])
    return campaign,attempt


def test_pilot_review_accepts_only_sampled_owned_resource_record(tmp_path,monkeypatch):
    campaign,attempt=fixture(tmp_path,monkeypatch)
    result=module.review(campaign,attempt,tmp_path/'review.json')
    assert result['status']=='NDI_neutral_UMA_pilot_independently_checked'
    assert result['peak_process_GPU_MiB']==50 and not result['experimental_accuracy_inferred']


@pytest.mark.parametrize('mutation',['foreign','missing_samples','wrong_gpu','model_error'])
def test_pilot_review_fails_closed(tmp_path,monkeypatch,mutation):
    campaign,attempt=fixture(tmp_path,monkeypatch)
    if mutation=='foreign': (attempt/'residency_guard/state.json').write_text(json.dumps({'foreign_resident_PIDs':[5],'own_residency_observed':True}))
    if mutation in ('missing_samples','wrong_gpu'):
        p=attempt/'gpu_memory/status.json';r=json.loads(p.read_text());r['own_samples']=0 if mutation=='missing_samples' else 2;r['GPU']='9' if mutation=='wrong_gpu' else '1';p.write_text(json.dumps(r))
    if mutation=='model_error': monkeypatch.setattr(module.uma,'cache_errors',lambda *a:['unconverged'])
    with pytest.raises(ValueError): module.review(campaign,attempt,tmp_path/'review.json')
