"""Receipt provenance checks reject swapped identities and corrupted scalars."""
import json
from pathlib import Path
import pytest

from redox.core.hashing import sha256_file
from redox.validation import sampling_receipts as receipts


def setup(tmp_path,monkeypatch):
    case=tmp_path/'case';case.mkdir();raw=case/'raw';raw.mkdir();namespace=tmp_path/'namespace.json'
    namespace.write_text('{}');task=dict(key='test',charge=0,mult=1)
    protocol=dict(xc='test-only');model=dict(name='test-only')
    plan=dict(tasks=[task],raw_root=str(raw),experimental_model=model,electronic_protocol=protocol)
    (case/'plan.json').write_text(json.dumps(plan))
    energy=dict(energy_protocol=protocol,experimental_method_id='test-model',converged_smd=True,converged_gas=True,
        optimized=False,charge=0,mult=1,e_smd_eV=-20.,e_gas_eV=-19.)
    (raw/'experimental_energy.json').write_text(json.dumps(energy))
    report=dict(status='corrected_sampling_energy_independently_checked',key='test',namespace=str(namespace),
        portable_method_id='test-model',actual_checkpoint_frames_verified=True,occupied_wavefunctions_checked=True,
        native_G_qualified=False,primary_candidate_admission=False,task=task,experimental_model=model,
        electronic_protocol=protocol,E_smd_eV=-20.,E_gas_eV=-19.,source_sha256={str(p):sha256_file(p)
            for p in (namespace,case/'plan.json',raw/'experimental_energy.json')})
    path=tmp_path/'review.json';path.write_text(json.dumps(report))
    monkeypatch.setattr(receipts,'verify',lambda p:dict(case=str(case)))
    return path,namespace,report


def test_cached_values_keep_their_exact_raw_record_and_context(tmp_path,monkeypatch):
    path,namespace,_=setup(tmp_path,monkeypatch)
    assert receipts.checked_energy(path,namespace,'test','test-model')['E_smd_eV']==-20.


@pytest.mark.parametrize('bad',['scalar','key','namespace','model','protocol','missing_raw_binding','changed_raw'])
def test_cached_receipt_cannot_silently_change_a_property(tmp_path,monkeypatch,bad):
    path,namespace,data=setup(tmp_path,monkeypatch)
    if bad=='scalar':data['E_smd_eV']=-10.
    if bad=='key':data['key']='another'
    if bad=='namespace':data['namespace']='/different/namespace.json'
    if bad=='model':data['experimental_model']={'name':'different'}
    if bad=='protocol':data['electronic_protocol']={'xc':'different'}
    if bad=='missing_raw_binding':data['source_sha256'].pop(str(tmp_path/'case/raw/experimental_energy.json'))
    if bad=='changed_raw':(tmp_path/'case/raw/experimental_energy.json').write_text('{}')
    path.write_text(json.dumps(data))
    with pytest.raises(ValueError):receipts.checked_energy(path,namespace,'test','test-model')
