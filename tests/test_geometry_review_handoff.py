import importlib.util
import json
from pathlib import Path
import pytest

ROOT=Path(__file__).resolve().parents[1]
s=importlib.util.spec_from_file_location('geometry_handoff',ROOT/'scripts/pipeline/publish_geometry_review_handoff.py');m=importlib.util.module_from_spec(s);s.loader.exec_module(m)


def fixture(tmp_path,monkeypatch):
    monkeypatch.setattr(m,'ROOT',tmp_path)
    source=tmp_path/'source.json';target=tmp_path/'results/review.json';campaign=tmp_path/'case'
    data=dict(status='alternate_smooth_geometry_and_resource_pilot_independently_checked',campaign=str(campaign),primary_candidate_admission=False,active_electronic_energy_qualified=False,native_Hessian_or_G_qualified=False,source_sha256={})
    source.write_text(json.dumps(data));return source,target,campaign,data


def test_exact_review_handoff_is_idempotent_and_immutable(tmp_path,monkeypatch):
    source,target,campaign,_=fixture(tmp_path,monkeypatch)
    assert m.validated_copy(source,target,campaign)['status']=='published'
    assert target.read_bytes()==source.read_bytes()
    assert m.validated_copy(source,target,campaign)['status']=='already_identical'
    target.write_text('other evidence')
    with pytest.raises(ValueError,match='different existing'):m.validated_copy(source,target,campaign)
    assert target.read_text()=='other evidence'


@pytest.mark.parametrize('field,value',[('campaign','/other'),('native_Hessian_or_G_qualified',True),('active_electronic_energy_qualified',True)])
def test_wrong_identity_or_promoted_scope_cannot_be_published(tmp_path,monkeypatch,field,value):
    source,target,campaign,data=fixture(tmp_path,monkeypatch);data[field]=value;source.write_text(json.dumps(data))
    with pytest.raises(ValueError):m.validated_copy(source,target,campaign)
    assert not target.exists()
