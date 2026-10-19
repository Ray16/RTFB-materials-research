"""Both readers must reject unqualified or source-changed PCM memory overrides."""
import hashlib,json
import pytest
from redox.qm import hessian_memory as M,native_thermal as N
from test_independent_native_audit import audit_module


def fixture(tmp_path,monkeypatch):
    source=tmp_path/'proof.json';source.write_text('immutable full-matrix PCM proof')
    data=dict(status='execution_pcm_curvature_qualified',protocol=N.NATIVE_THERMAL,
        software={'pyscf':'2.14.0','gpu4pyscf':'1.8.1'},mode=M.PCM_MODE,
        scientific_model_changed=False,surface_grid_changed=False,adapter_sha256='adapter',
        cpu_pcm_library_sha256='cpu',gpu_pcm_library_sha256='gpu',
        isolated_proof_relative_path='proof.json',isolated_proof_sha256=M.digest(source),
        controls=[dict(mult=m,native_gates_passed=True,relative_matrix_error=1e-10,
           maximum_matrix_difference=1e-10,maximum_thermal_difference_eV=1e-8,
           evidence_sha256={'proof.json':M.digest(source)}) for m in [1,2]])
    data['content_sha256']=hashlib.sha256(json.dumps(data,sort_keys=True,separators=(',',':')).encode()).hexdigest()
    cert=tmp_path/'config/native_pcm_curvature_memory_validation.json';cert.parent.mkdir();cert.write_text(json.dumps(data));monkeypatch.setattr(M,'ROOT',tmp_path);monkeypatch.setattr(M,'PCM_CERTIFICATE',cert)
    s=dict(mode=M.PCM_MODE,qualification_only=False,production_admission=True,
       scientific_model_changed=False,surface_grid_changed=False,
       qualification_sha256=data['content_sha256'],source_geometry_sha256='geom',
       adapter_sha256='adapter',cpu_pcm_library_sha256='cpu',gpu_pcm_library_sha256='gpu',
       isolated_proof_sha256=M.digest(source),target_evidence_sha256={'proof.json':M.digest(source)})
    r=dict(identity={'geometry_sha256':'geom'},software=data['software'],execution_pcm_curvature=s)
    audit=audit_module();monkeypatch.setattr(audit,'ROOT',tmp_path)
    return source,r,lambda x:M.record_pcm_ok(x,N.NATIVE_THERMAL),audit.native_execution_pcm_ok


@pytest.mark.parametrize('change',['qualification','model','surface','adapter','library','source','seal','missing_target_evidence','changed_evidence'])
def test_both_pcm_readers_reject_wrong_execution_or_stale_evidence(tmp_path,monkeypatch,change):
    source,r,production,independent=fixture(tmp_path,monkeypatch)
    assert production(r) and independent(r) and production({}) and independent({})
    s=r['execution_pcm_curvature']
    if change=='qualification':s['qualification_only']=True
    if change=='model':s['scientific_model_changed']=True
    if change=='surface':s['surface_grid_changed']=True
    if change=='adapter':s['adapter_sha256']='different'
    if change=='library':s['cpu_pcm_library_sha256']='different'
    if change=='source':s['source_geometry_sha256']='different'
    if change=='seal':s['qualification_sha256']='different'
    if change=='missing_target_evidence':s['target_evidence_sha256']={}
    if change=='changed_evidence':source.write_text('changed numerical proof')
    assert not production(r) and not independent(r)
