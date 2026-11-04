import json
from pathlib import Path
import pytest
from redox.workflows import charged_native_followup as followup

FILES=('worker.py','helper.py','kernel.cu','checkpoints.py','response_helper.py','restart_helper.py','response_resume_helper.py')


def setup(tmp_path,monkeypatch):
    specs={}
    for name,nao in [('pilot',1100),('child',881)]:
        root=tmp_path/name;case=root/'case';sources=root/'sources';case.mkdir(parents=True);sources.mkdir()
        for leaf in FILES:(sources/leaf).write_text('test-only shared source')
        plan=dict(charged_polymer_native_component=True,response_nao=nao,response_nocc=177,mult=1,charge=0 if name=='child' else 2,solvent_model={'name':'test'},numerical_settings={'step':.001})
        (case/'plan.json').write_text(json.dumps(plan));(case/'thermodynamic_component.json').write_text(json.dumps(dict(G_by_cutoff_eV={c:-100. for c in ('25','50','100')})))
        ns=root/'namespace.json';ns.write_text('{}');specs[str(ns)]=dict(case=str(case),sources=str(sources))
    monkeypatch.setattr(followup,'verify',lambda p:specs[str(p)])
    ns=tmp_path/'pilot/namespace.json';review=dict(status='charged_polymer_native_component_independently_reconstructed',native_G_qualified=True,namespace=str(ns),source_sha256={},charge=2,mult=1,G_by_cutoff_eV={c:-100. for c in ('25','50','100')},peak_GPU_MiB=23000,maximum_column_error=1e-5,minimum_frequency_cm=8.)
    path=tmp_path/'review.json';path.write_text(json.dumps(review))
    manifest=dict(prerequisite_namespace=str(ns),prepared=[dict(key='child',namespace=str(tmp_path/'child/namespace.json'))])
    return manifest,path,review


def test_completed_matching_pilot_releases_only_its_declared_class(tmp_path,monkeypatch):
    manifest,path,_=setup(tmp_path,monkeypatch)
    assert len(followup.prerequisite(manifest,path))==1


@pytest.mark.parametrize('bad',['pending','incomplete_G','larger_dimension','different_kernel'])
def test_pending_or_incompatible_pilot_cannot_release_followups(tmp_path,monkeypatch,bad):
    manifest,path,review=setup(tmp_path,monkeypatch)
    if bad=='pending':review['native_G_qualified']=False
    if bad=='incomplete_G':review['G_by_cutoff_eV'].pop('100')
    if bad=='larger_dimension':
        p=tmp_path/'child/case/plan.json';x=json.loads(p.read_text());x['response_nao']=1300;p.write_text(json.dumps(x))
    if bad=='different_kernel':(tmp_path/'child/sources/kernel.cu').write_text('changed')
    path.write_text(json.dumps(review))
    with pytest.raises(ValueError):followup.prerequisite(manifest,path)


def larger_setup(tmp_path,monkeypatch):
    from redox.workflows.native_component_resources import LARGE_CLASS,larger_profile
    manifest,path,review=setup(tmp_path,monkeypatch)
    for name in ('pilot','child'):
        p=tmp_path/name/'case/plan.json';data=json.loads(p.read_text())
        data.update(atoms=137,response_nao=1319,response_nocc=249,charge=4,
                    prepared_case={'occupied_sites':[2,3]},maximum_host_RSS_GiB=80,
                    maximum_wall_seconds=28800,minimum_initial_storage_GiB=38,minimum_live_storage_GiB=32,
                    **larger_profile(137,1319,249,4,1,[2,3]))
        p.write_text(json.dumps(data));(tmp_path/name/'sources/native_resource_helper.py').write_text('same frozen resource bounds')
    manifest['native_resource_class']=LARGE_CLASS
    review.update(charge=4,native_resource_class=LARGE_CLASS,peak_GPU_MiB=45000)
    path.write_text(json.dumps(review));return manifest,path,review


def test_larger_class_requires_and_accepts_its_own_passed_pilot(tmp_path,monkeypatch):
    manifest,path,_=larger_setup(tmp_path,monkeypatch)
    assert len(followup.prerequisite(manifest,path))==1


@pytest.mark.parametrize('bad',['small_pilot','different_helper','lower_memory_gate'])
def test_larger_followups_reject_unproven_pilot_or_changed_resource_class(tmp_path,monkeypatch,bad):
    manifest,path,review=larger_setup(tmp_path,monkeypatch)
    if bad=='small_pilot':review.update(charge=2,native_resource_class='original_32GiB_class')
    if bad=='different_helper':(tmp_path/'child/sources/native_resource_helper.py').write_text('changed limits')
    if bad=='lower_memory_gate':
        p=tmp_path/'child/case/plan.json';data=json.loads(p.read_text());data['minimum_GPU_memory_MiB']=31744;p.write_text(json.dumps(data))
    path.write_text(json.dumps(review))
    with pytest.raises(ValueError):followup.prerequisite(manifest,path)
