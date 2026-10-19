import hashlib
import json
import pytest
from redox.qm import hessian_memory as M
from redox.qm.native_thermal import NATIVE_THERMAL


def make_certificate(tmp_path, monkeypatch):
    source = tmp_path/'control.json'
    source.write_text('independent control evidence')
    controls = [dict(mult=m, relative_matrix_error=1e-9, maximum_matrix_difference=1e-10,
                     maximum_thermal_difference_eV=1e-10, native_gates_passed=True,
                     evidence_sha256={'control.json':M.digest(source)}) for m in (1,2)]
    data = dict(status='execution_budget_qualified',protocol=NATIVE_THERMAL,
                software={'pyscf':'2.14.0','gpu4pyscf':'1.8.1'},caps_GiB=[24],
                scientific_model_changed=False,controls=controls,allocation_policy=M.ALLOCATION_POLICY)
    data['content_sha256'] = hashlib.sha256(json.dumps(data,sort_keys=True,separators=(',',':')).encode()).hexdigest()
    path = tmp_path/'certificate.json'
    path.write_text(json.dumps(data))
    monkeypatch.setattr(M,'ROOT',tmp_path)
    monkeypatch.setattr(M,'CERTIFICATE',path)
    return source, path


def test_larger_budget_fails_closed_without_controls(tmp_path,monkeypatch):
    monkeypatch.setattr(M,'CERTIFICATE',tmp_path/'absent.json')
    assert M.record_budget_ok({},NATIVE_THERMAL)
    with pytest.raises(OSError):M.execution_budget(NATIVE_THERMAL,24)
    assert not M.record_budget_ok({'execution_memory':dict(cap_GiB=24)},NATIVE_THERMAL)


def test_qualified_override_preserves_the_model_and_rejects_changed_evidence(tmp_path,monkeypatch):
    source,path = make_certificate(tmp_path,monkeypatch)
    budget = M.execution_budget(NATIVE_THERMAL,24)
    assert budget['cap_GiB']==24 and budget['scientific_model_changed'] is False
    assert M.record_budget_ok({'execution_memory':budget},NATIVE_THERMAL)
    source.write_text('changed evidence')
    assert not M.record_budget_ok({'execution_memory':budget},NATIVE_THERMAL)


def test_qualification_cannot_admit_production_and_unsupported_budgets_fail():
    budget = M.execution_budget(NATIVE_THERMAL,24,qualification=True)
    assert not M.record_budget_ok({'execution_memory':budget},NATIVE_THERMAL)
    with pytest.raises(ValueError):M.execution_budget(NATIVE_THERMAL,96)


def test_default_budget_does_not_require_a_new_certificate(monkeypatch):
    monkeypatch.delenv('REDOX_HESSIAN_MEMORY_GIB',raising=False)
    budget=M.execution_budget(NATIVE_THERMAL)
    assert budget['cap_GiB']==6
    assert M.record_budget_ok({'execution_memory':budget},NATIVE_THERMAL)


def test_larger_response_workspace_does_not_inflate_integral_batches():
    budget=M.execution_budget(NATIVE_THERMAL,24,qualification=True)
    assert M.stage_budget_bytes(budget,'gpu4pyscf.df.hessian.uhf','_get_veff')==24*1024**3
    assert M.stage_budget_bytes(budget,'gpu4pyscf.df.hessian.uhf','_jk_energy_per_atom')==6*1024**3
    assert M.stage_budget_bytes(budget,'gpu4pyscf.df.hessian.rhf','_j_energy_per_atom')==6*1024**3


def test_distinct_sixteen_budget_certificate_preserves_twentyfour(tmp_path,monkeypatch):
    source,path24=make_certificate(tmp_path,monkeypatch)
    original24=path24.read_bytes();data=json.loads(path24.read_text());data.pop('content_sha256')
    data.update(caps_GiB=[16],allocation_policy=M.ALLOCATION_POLICY_16)
    data['content_sha256']=hashlib.sha256(json.dumps(data,sort_keys=True,separators=(',',':')).encode()).hexdigest()
    path16=tmp_path/'config/native_hessian_memory_16_validation.json';path16.parent.mkdir();path16.write_text(json.dumps(data));monkeypatch.setattr(M,'CERTIFICATE_16',path16)
    budget=M.execution_budget(NATIVE_THERMAL,16)
    assert M.record_budget_ok({'execution_memory':budget},NATIVE_THERMAL)
    assert M.stage_budget_bytes(budget,'gpu4pyscf.df.hessian.uhf','_get_veff')==16*1024**3
    assert M.stage_budget_bytes(budget,'gpu4pyscf.df.hessian.uhf','_jk_energy_per_atom')==6*1024**3
    assert path24.read_bytes()==original24 and M.execution_budget(NATIVE_THERMAL,24)['cap_GiB']==24
    from test_independent_native_audit import audit_module
    audit=audit_module();monkeypatch.setattr(audit,'ROOT',tmp_path)
    record={'execution_memory':budget,'software':{'pyscf':'2.14.0','gpu4pyscf':'1.8.1'}}
    assert audit.native_execution_budget_ok(record)
    source.write_text('changed evidence')
    assert not M.record_budget_ok(record,NATIVE_THERMAL)
    assert not audit.native_execution_budget_ok(record)
