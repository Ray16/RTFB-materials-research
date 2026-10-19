"""Both readers must reject unqualified CDS variants and stale contact evidence."""
import copy,hashlib,json
import pytest
from redox.qm import hessian_memory as M,native_thermal as N
from test_independent_native_audit import audit_module,fixture_record


def proof(tmp_path,monkeypatch,protocol,step=1e-5):
    source=tmp_path/'evidence.json';source.write_text('immutable replay and target evidence')
    data=dict(schema=1,status='execution_derivatives_qualified',protocol=protocol,
              software={'pyscf':'2.14.0','gpu4pyscf':'1.8.1'},cds_step_Bohr=step,
              library_default_cds_step_Bohr=1e-4,scientific_model_changed=False,
              worst_contact_validation=dict(steps_Bohr=[.0005,.001],relative_acceptance=.01,symmetry_acceptance=1e-6),
              controls=[dict(mult=m,native_gates_passed=True,relative_matrix_error=1e-10,
                maximum_matrix_difference=1e-10,maximum_thermal_difference_eV=1e-9,
                evidence_sha256={'evidence.json':M.digest(source)}) for m in [1,2]])
    if step==1e-6:data.update(target_native_gates_passed=True,qualification_evidence_sha256={'evidence.json':M.digest(source)})
    data['content_sha256']=hashlib.sha256(json.dumps(data,sort_keys=True,separators=(',',':')).encode()).hexdigest()
    filename='native_cds_step_validation.json' if step==1e-5 else 'native_cds_step_1e6_validation.json'
    cert=tmp_path/'config'/filename;cert.parent.mkdir(exist_ok=True);cert.write_text(json.dumps(data))
    monkeypatch.setattr(M,'ROOT',tmp_path);monkeypatch.setattr(M,'CDS_CERTIFICATE',cert)
    monkeypatch.setattr(M,'CDS_CERTIFICATE_1E6',cert if step==1e-6 else tmp_path/'config/native_cds_step_1e6_validation.json')
    settings=dict(cds_step_Bohr=step,qualification_only=False,scientific_model_changed=False,
                  qualification_sha256=data['content_sha256'],source_geometry_sha256='geom',
                  worst_contact_indices=[0,1],worst_contact_checks=[dict(column=i,step_Bohr=s,relative_error=.0001) for i in [0,1] for s in [.0005,.001]],
                  target_evidence_sha256={'evidence.json':M.digest(source)})
    return source,settings


@pytest.mark.parametrize('step',[1e-5,1e-6])
def test_each_derivative_step_uses_its_own_certificate(tmp_path,monkeypatch,step):
    audit=audit_module();monkeypatch.setattr(audit,'ROOT',tmp_path)
    source,settings=proof(tmp_path,monkeypatch,N.NATIVE_THERMAL,step)
    record=dict(identity={'geometry_sha256':'geom'},software={'pyscf':'2.14.0','gpu4pyscf':'1.8.1'},execution_derivatives=settings)
    assert M.record_derivatives_ok(record,N.NATIVE_THERMAL,2) and audit.native_execution_derivatives_ok(record,2)
    if step==1e-6:
        (tmp_path/'config/native_cds_step_1e6_validation.json').unlink()
        assert not M.record_derivatives_ok(record,N.NATIVE_THERMAL,2)
        assert not audit.native_execution_derivatives_ok(record,2)


@pytest.mark.parametrize('change',['target_not_verified','missing_qualification_evidence'])
def test_finer_step_rejects_incomplete_target_qualification(tmp_path,monkeypatch,change):
    audit=audit_module();monkeypatch.setattr(audit,'ROOT',tmp_path)
    source,settings=proof(tmp_path,monkeypatch,N.NATIVE_THERMAL,1e-6)
    record=dict(identity={'geometry_sha256':'geom'},software={'pyscf':'2.14.0','gpu4pyscf':'1.8.1'},execution_derivatives=settings)
    path=tmp_path/'config/native_cds_step_1e6_validation.json';data=json.loads(path.read_text());data.pop('content_sha256')
    if change=='target_not_verified':data['target_native_gates_passed']=False
    else:data['qualification_evidence_sha256']={}
    data['content_sha256']=hashlib.sha256(json.dumps(data,sort_keys=True,separators=(',',':')).encode()).hexdigest();settings['qualification_sha256']=data['content_sha256'];path.write_text(json.dumps(data))
    assert not M.record_derivatives_ok(record,N.NATIVE_THERMAL,2)
    assert not audit.native_execution_derivatives_ok(record,2)


def readers(tmp_path,monkeypatch):
    audit=audit_module();monkeypatch.setattr(audit,'ROOT',tmp_path)
    source,settings=proof(tmp_path,monkeypatch,N.NATIVE_THERMAL)
    r=dict(identity={'geometry_sha256':'geom'},software={'pyscf':'2.14.0','gpu4pyscf':'1.8.1'},execution_derivatives=settings)
    return source,r,lambda x:M.record_derivatives_ok(x,N.NATIVE_THERMAL,2),lambda x:audit.native_execution_derivatives_ok(x,2)


@pytest.mark.parametrize('change',['missing_control_evidence','qualification_only','bad_seal','changed_model','wrong_step','wrong_geometry','missing_contact_step','large_contact_error','nonfinite_error','out_of_range_column','missing_target_evidence'])
def test_both_readers_reject_incomplete_or_stale_derivative_proof(tmp_path,monkeypatch,change):
    source,r,production,independent=readers(tmp_path,monkeypatch)
    assert production(r) and independent(r)
    s=r['execution_derivatives']
    if change=='missing_control_evidence':source.unlink()
    if change=='qualification_only':s['qualification_only']=True
    if change=='bad_seal':s['qualification_sha256']='wrong'
    if change=='changed_model':s['scientific_model_changed']=True
    if change=='wrong_step':s['cds_step_Bohr']=1e-6
    if change=='wrong_geometry':s['source_geometry_sha256']='different geometry'
    if change=='missing_contact_step':s['worst_contact_checks'].pop()
    if change=='large_contact_error':s['worst_contact_checks'][0]['relative_error']=.02
    if change=='nonfinite_error':s['worst_contact_checks'][0]['relative_error']=float('nan')
    if change=='out_of_range_column':s['worst_contact_indices']=[0,6]
    if change=='missing_target_evidence':s['target_evidence_sha256']={}
    assert not production(r) and not independent(r)


def test_both_readers_keep_historical_records_and_reject_changed_evidence(tmp_path,monkeypatch):
    source,r,production,independent=readers(tmp_path,monkeypatch)
    assert production({}) and independent({})
    source.write_text('changed target or control output')
    assert not production(r) and not independent(r)


def test_independent_full_admission_requires_contact_proof(tmp_path,monkeypatch):
    audit=audit_module();path,r,raw=fixture_record(tmp_path,audit)
    source,settings=proof(tmp_path,monkeypatch,audit.NATIVE_THERMAL_CONTRACT)
    monkeypatch.setattr(audit,'ROOT',tmp_path)
    settings['source_geometry_sha256']=r['identity']['geometry_sha256'];r['execution_derivatives']=settings;path.write_text(json.dumps(r))
    assert audit._surface_thermal_inputs(tmp_path,raw)[0]==1.
    settings['worst_contact_checks']=[];path.write_text(json.dumps(r))
    assert audit._surface_thermal_inputs(tmp_path,raw)[0] is None


def test_production_full_admission_requires_contact_proof(tmp_path,monkeypatch):
    from test_native_thermal import record
    monkeypatch.setattr(N,'implementation_controls',lambda:[{'verified':'test controls'}])
    (tmp_path/'opt.xyz').write_text('3\n\nO 0 0 0\nH .7 0 .5\nH -.7 0 .5\n')
    path,r=record(tmp_path);source,settings=proof(tmp_path,monkeypatch,N.NATIVE_THERMAL)
    settings['source_geometry_sha256']=r['identity']['geometry_sha256'];r['execution_derivatives']=settings;path.write_text(json.dumps(r))
    assert N._load_native_surface_record(tmp_path,-1,2)
    settings['worst_contact_checks']=[];path.write_text(json.dumps(r))
    assert N._load_native_surface_record(tmp_path,-1,2) is None
