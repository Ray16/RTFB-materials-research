import json
from pathlib import Path
import pytest
from redox.core.hashing import sha256_file
from redox.workflows import relocate_unstarted_namespace as relocation


def setup(tmp_path,monkeypatch):
    old=tmp_path/'old';case=old/'case';sources=old/'sources';case.mkdir(parents=True);sources.mkdir()
    data=sources/'input.xyz';data.write_text('test-only input')
    dependency=tmp_path/'node_local.py';dependency.write_text('identical dependency')
    shared=tmp_path/'shared.py';shared.write_bytes(dependency.read_bytes())
    plan=dict(source_id='test',state='rest',source_sha256={str(data):sha256_file(data),str(dependency):sha256_file(dependency)},input=str(data),charge=2)
    (case/'plan.json').write_text(json.dumps(plan));spec=dict(case=str(case),sources=str(sources),archive=str(tmp_path/'archive'),plan_sha256=sha256_file(case/'plan.json'))
    namespace=old/'namespace.json';namespace.write_text(json.dumps(spec))
    def verify(path):
        current=json.loads(path.read_text());p=Path(current['case'])/'plan.json'
        assert sha256_file(p)==current['plan_sha256']
        for name,wanted in json.loads(p.read_text())['source_sha256'].items():assert sha256_file(name)==wanted
        return current
    monkeypatch.setattr(relocation,'verify',verify)
    return namespace,dependency,shared


def test_only_locations_change_and_original_is_preserved(tmp_path,monkeypatch):
    namespace,dependency,shared=setup(tmp_path,monkeypatch);before=namespace.read_bytes()
    result=relocation.relocate(namespace,tmp_path/'new',{str(dependency):str(shared)})
    assert namespace.read_bytes()==before and result['scientific_parameters_changed'] is False
    plan=json.loads((tmp_path/'new/case/plan.json').read_text())
    assert plan['charge']==2 and str(shared) in plan['source_sha256']
    assert Path(plan['input']).read_text()=='test-only input'


@pytest.mark.parametrize('bad',['changed_bytes','already_started'])
def test_changed_dependencies_or_started_calculations_are_rejected(tmp_path,monkeypatch,bad):
    namespace,dependency,shared=setup(tmp_path,monkeypatch)
    if bad=='changed_bytes':shared.write_text('different')
    else:(namespace.parent/'case/status.json').write_text('{}')
    with pytest.raises(ValueError):relocation.relocate(namespace,tmp_path/'new',{str(dependency):str(shared)})
