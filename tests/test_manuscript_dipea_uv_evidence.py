import hashlib
import importlib.util
import json
from pathlib import Path
import shutil

import pytest

ROOT=Path(__file__).parents[1]
spec=importlib.util.spec_from_file_location('checked_manuscript_dipea_uv',ROOT/'manuscript/check_claims.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)


@pytest.fixture
def relocated_epoch(tmp_path,monkeypatch):
    base=tmp_path/'manuscript/evidence';base.mkdir(parents=True)
    source=ROOT/'manuscript/evidence/2026-10-26_dipea_uv';target=base/'review'
    shutil.copytree(source,target)
    (base/'current.json').write_text(json.dumps({'manifest':'manuscript/evidence/review/manifest.json'}))
    monkeypatch.setattr(module,'__file__',str(tmp_path/'manuscript/check_claims.py'))
    return target


def test_relocated_claims_reconstruct_uv_without_promoting_hybrid(relocated_epoch):
    result=module.audit()['dipea_uv_anchor']
    assert result['source_logKH']==18.64
    assert result['capture_difference_kJ_mol']==pytest.approx(-23.37063570)
    assert result['relative_difference_kJ_mol']==pytest.approx(-.56283794)
    assert not result['hybrid_is_calibrated_prediction']


@pytest.mark.parametrize('alteration',['source_value','hybrid_value','promoted_prediction','unmatched_source_hash'])
def test_uv_claim_mismatch_or_promotion_is_rejected(relocated_epoch,alteration):
    path=relocated_epoch/'dipea_uv_anchor.json';row=json.loads(path.read_text())
    if alteration=='source_value':row['source']['reported_value']=18.60
    elif alteration=='hybrid_value':row['candidate_hybrid_sensitivities'][0]['hybrid_capture_sensitivity_kJ_mol']+=1
    elif alteration=='promoted_prediction':row['candidate_hybrid_sensitivities'][0]['calibrated_reaction_prediction']=True
    else:row['source']['technique']='altered'
    path.write_text(json.dumps(row))
    if alteration!='unmatched_source_hash':
        review_path=relocated_epoch/'dipea_uv_review.json';review=json.loads(review_path.read_text())
        key=next(k for k in review['source_sha256'] if Path(k).name=='summary.json')
        review['source_sha256'][key]=hashlib.sha256(path.read_bytes()).hexdigest();review_path.write_text(json.dumps(review))
    with pytest.raises(ValueError,match='UV|uv'):module.audit()
