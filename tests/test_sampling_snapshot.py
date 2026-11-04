"""Snapshot publication preserves exact parents and never fills pending slots."""
from copy import deepcopy
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from redox.core.hashing import sha256_file
from redox.validation import sampling_snapshot as api


def fixture(tmp_path,monkeypatch):
    keys=['all_same__occ_2__conf'+str(i).zfill(3) for i in range(1,241)]
    manifest=tmp_path/'manifest.json'
    manifest.write_text(json.dumps(dict(target_states=240,experimental_method_id='method',states=[dict(key=k) for k in keys])))
    case=tmp_path/'case';case.mkdir();geometry=case/'input.xyz';geometry.write_text('exact geometry')
    ns=tmp_path/'namespace.json';ns.write_text(json.dumps(dict(case=str(case))))
    ereview=tmp_path/'energy_review.json';ereview.write_text('checked electronic fixture')
    nreview=tmp_path/'native_review.json';nreview.write_text('checked native fixture')
    task=dict(key=keys[0],charge=0,mult=1)
    e=dict(task=task,namespace=str(ns),experimental_model=dict(name='model'),electronic_protocol=dict(xc='method'),
        E_smd_eV=-2.,E_gas_eV=-1.9)
    n=dict(task=task,experimental_model=e['experimental_model'],electronic_protocol=e['electronic_protocol'],
        native_settings=dict(xc='native'),thermal_protocol=dict(pressure='declared'),standard_state='declared',
        temperature_K=298.15,charge=0,multiplicity=1,geometry_sha256=sha256_file(geometry),experimental_method_id='method',
        electronic_E_smd_eV=-2.,nuclear_by_cutoff_eV={k:.1 for k in ('25','50','100')},
        G_by_cutoff_eV={k:-1.9 for k in ('25','50','100')},minimum_frequency_cm=10.,maximum_column_error=.0001,
        observed_peak_GPU_MiB=20000,actual_SCF_frames_checked=17)
    erows=[dict(key=k,status='pending',G_eV=None,native_G_qualified=False,E_smd_eV=None,E_gas_eV=None) for k in keys]
    erows[0].update(status='electronic_independently_checked',namespace=str(ns),review=str(ereview),
        review_sha256=sha256_file(ereview),E_smd_eV=e['E_smd_eV'],E_gas_eV=e['E_gas_eV'])
    nrows=[dict(key=k,status='native_not_prepared',native_G_qualified=False,G_by_cutoff_eV=None,nuclear_by_cutoff_eV=None) for k in keys]
    nrows[0].update(status='native_independently_checked',native_G_qualified=True,namespace=str(ns),review=str(nreview),
        review_sha256=sha256_file(nreview),G_by_cutoff_eV=deepcopy(n['G_by_cutoff_eV']),nuclear_by_cutoff_eV=deepcopy(n['nuclear_by_cutoff_eV']),
        electronic_E_smd_eV=n['electronic_E_smd_eV'],context=deepcopy({k:n[k] for k in api.CONTEXT}))
    base=dict(declared_inputs=240,experimental_method_id='method',independent_review_errors=[],
        occupied_ensemble_cycles_qualified=False,pending_values_are_null=True,source_sha256={str(manifest):sha256_file(manifest)})
    energy=dict(base,status='complete_roster_corrected_electronic_review_snapshot',states=erows,checked_electronic_inputs=1,native_G_qualified=False)
    native=dict(base,status='complete_roster_corrected_native_review_snapshot',states=nrows,checked_native_inputs=1)
    ep=tmp_path/'electronic.json';np=tmp_path/'native.json'
    ep.write_text(json.dumps(energy));np.write_text(json.dumps(native))
    monkeypatch.setattr(api,'checked_energy',lambda *a:deepcopy(e))
    monkeypatch.setattr(api,'ReceiptCache',lambda:SimpleNamespace(read=lambda *a:deepcopy(n)))
    monkeypatch.setattr(api,'verify',lambda p:dict(case=str(case)))
    return ep,np,energy,native,n,e


def test_only_checked_components_are_published_and_pending_slots_stay_visible(tmp_path,monkeypatch):
    ep,np,*_=fixture(tmp_path,monkeypatch);report=api.audit(ep,np)
    assert report['checked_electronic_inputs']==report['checked_native_inputs']==1
    assert report['actual_SCF_frames_represented']==17
    assert not report['original_physical_QC_repeated'] and not report['occupied_ensemble_cycles_qualified']


@pytest.mark.parametrize('bad',['pending_energy','pending_G','native_value','native_context','census','missing_parent','missing_slot','source_roster','source_change'])
def test_snapshot_count_context_and_pending_values_cannot_be_promoted(tmp_path,monkeypatch,bad):
    ep,np,e,n,_,_=fixture(tmp_path,monkeypatch)
    if bad=='pending_energy':e['states'][1]['E_smd_eV']=0.
    if bad=='pending_G':n['states'][1]['G_by_cutoff_eV']={'50':0.}
    if bad=='native_value':n['states'][0]['G_by_cutoff_eV']['50']=-2.1
    if bad=='native_context':n['states'][0]['context']['temperature_K']=300
    if bad=='census':n['checked_native_inputs']=2
    if bad=='missing_parent':e['states'][0].update(status='pending',E_smd_eV=None,E_gas_eV=None)
    if bad=='missing_slot':n['states'].pop()
    if bad=='source_roster':e['states'][1]['key']=n['states'][1]['key']='invented input'
    if bad=='source_change':Path(next(iter(e['source_sha256']))).write_text('{}')
    ep.write_text(json.dumps(e));np.write_text(json.dumps(n))
    with pytest.raises(ValueError):api.audit(ep,np)


def test_separate_repeat_is_compared_without_replacing_the_native_parent(tmp_path,monkeypatch):
    _,_,_,_,n,e=fixture(tmp_path,monkeypatch);e['E_smd_eV']+=5.8e-11
    agreement=api.parent_agreement(n,e,n['geometry_sha256'])
    assert agreement['absolute_energy_difference_eV']>0 and not agreement['native_electronic_parent_replaced']
    for bad in ('energy','geometry','model','task'):
        row=deepcopy(e);geometry=n['geometry_sha256']
        if bad=='energy':row['E_smd_eV']+=1e-4
        if bad=='geometry':geometry='another geometry'
        if bad=='model':row['experimental_model']={}
        if bad=='task':row['task']['charge']=1
        with pytest.raises(ValueError):api.parent_agreement(n,row,geometry)
