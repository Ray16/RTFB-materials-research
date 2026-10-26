import importlib.util
import json
from pathlib import Path
import shutil

import pytest

ROOT=Path(__file__).parents[1]
spec=importlib.util.spec_from_file_location('checked_manuscript_dipea',ROOT/'manuscript/check_claims.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)


@pytest.fixture
def relocated_epoch(tmp_path,monkeypatch):
    base=tmp_path/'manuscript/evidence';base.mkdir(parents=True)
    src=ROOT/'manuscript/evidence/2026-10-26_dipea_nmr';target=base/'review'
    shutil.copytree(src,target)
    (base/'current.json').write_text(json.dumps({'manifest':'manuscript/evidence/review/manifest.json'}))
    monkeypatch.setattr(module,'__file__',str(tmp_path/'manuscript/check_claims.py'))
    return target


def test_relocated_manuscript_claims_reconstruct_dipea_without_promoting_hybrid(relocated_epoch):
    result=module.audit()['dipea_experimental_anchor']
    assert result['source_pKa']==18.7
    assert result['capture_difference_kJ_mol']==pytest.approx(-23.02815513)
    assert result['relative_difference_kJ_mol']==pytest.approx(-.22035737)
    assert not result['hybrid_is_calibrated_prediction']


@pytest.mark.parametrize('alteration', ['source_pka','hybrid_value','promoted_prediction'])
def test_altered_experimental_claim_or_hybrid_promotion_is_rejected(relocated_epoch,alteration):
    path=relocated_epoch/'dipea_experimental_anchor.json';row=json.loads(path.read_text())
    if alteration=='source_pka':row['source']['reported_pKa']=18.6
    elif alteration=='hybrid_value':row['candidate_hybrid_sensitivities'][0]['hybrid_capture_sensitivity_kJ_mol']+=1
    else:row['candidate_hybrid_sensitivities'][0]['hybrid_is_calibrated_reaction_prediction']=True
    path.write_text(json.dumps(row))
    with pytest.raises(ValueError,match='DIPEA|Hybrid'):module.audit()
