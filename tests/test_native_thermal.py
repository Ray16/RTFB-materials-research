import json
import pytest
from redox.qm import native_thermal as N
from redox.qm.thermal import spin_degeneracy_eV


@pytest.fixture(autouse=True)
def controls(monkeypatch):
    monkeypatch.setattr(N,"implementation_controls",lambda: [{"verified":"test controls"}])


def record(sd,q=-1,m=2):
    p,ident=N.cache_path(sd,q,m);p.parent.mkdir(parents=True,exist_ok=True)
    spin=spin_degeneracy_eV(m)
    return p,dict(status='derivative_checked_solution_thermal',identity=ident,
        stationary_minimum=True,converged=True,grid_response=True,
        g_nuclear_qrrho_eV=2.,g_spin_partition_eV=spin,g_thermal_total_eV=2.+spin,
        nuclear_by_cutoff_eV={'50':2.},frequencies_cm=[50.,100.,200.],gradient_max=.0001,
        maximum_asymmetry=1e-10,
        implementation_controls=[{"verified":"test controls"}],
        derivative_checks=[dict(atom_index=i,axis=a,step_Bohr=s,relative_error=.0001)
            for i,a in ((0,1),(0,2)) for s in (.0005,.001)],
        smd_spin_diagnostics=dict(spin_contam=.01,homo_eV=-1.,s_squared=.76),software={'pyscf':'2.14.0','gpu4pyscf':'1.8.1'})


def test_valid_native_record_has_one_explicit_thermal_definition(tmp_path):
    (tmp_path/'opt.xyz').write_text('3\n\nO 0 0 0\nH .7 0 .5\nH -.7 0 .5\n')
    p,r=record(tmp_path);p.write_text(json.dumps(r))
    assert N._load_native_surface_record(tmp_path,-1,2)
    r['g_thermal_total_eV']=2.;p.write_text(json.dumps(r))
    assert N._load_native_surface_record(tmp_path,-1,2) is None


@pytest.mark.parametrize('field,value', [('gradient_max',float('nan')),
    ('grid_response',False),('frequencies_cm',[-12.,50.]),('converged',False),
    ('derivative_checks',[dict(relative_error=float('nan')),dict(relative_error=0.)]),
    ('smd_spin_diagnostics',dict(spin_contam=.5))])
def test_bad_native_physics_never_loads(tmp_path,field,value):
    (tmp_path/'opt.xyz').write_text('3\n\nO 0 0 0\nH .7 0 .5\nH -.7 0 .5\n')
    p,r=record(tmp_path);r[field]=value;p.write_text(json.dumps(r))
    assert N._load_native_surface_record(tmp_path,-1,2) is None


def test_native_cache_rejects_geometry_and_spin_changes(tmp_path):
    (tmp_path/'opt.xyz').write_text('3\n\nO 0 0 0\nH .7 0 .5\nH -.7 0 .5\n')
    p,r=record(tmp_path);p.write_text(json.dumps(r))
    assert N._load_native_surface_record(tmp_path,-1,2)
    assert N._load_native_surface_record(tmp_path,-1,4) is None
    (tmp_path/'opt.xyz').write_text('3\n\nO 0 0 0\nH .8 0 .5\nH -.7 0 .5\n')
    assert N._load_native_surface_record(tmp_path,-1,2) is None


@pytest.mark.parametrize("change", ["negative_gradient","negative_error","bad_asymmetry","mode_count","missing_steps","control_identity"])
def test_complete_numerical_contract(tmp_path,change):
    (tmp_path/'opt.xyz').write_text('3\n\nO 0 0 0\nH .7 0 .5\nH -.7 0 .5\n')
    p,r=record(tmp_path)
    if change=='negative_gradient':r['gradient_max']=-1.
    if change=='negative_error':r['derivative_checks'][0]['relative_error']=-1.
    if change=='bad_asymmetry':r['maximum_asymmetry']=.1
    if change=='mode_count':r['frequencies_cm']=[50.]
    if change=='missing_steps':r['derivative_checks']=r['derivative_checks'][:2]
    if change=='control_identity':r['implementation_controls']=[]
    p.write_text(json.dumps(r));assert N._load_native_surface_record(tmp_path,-1,2) is None


def test_policy_missing_native_never_uses_legacy(tmp_path):
    from redox.core.common import apply_thermal_policy,free_energy
    (tmp_path/'opt.xyz').write_text('3\n\nO 0 0 0\nH .7 0 .5\nH -.7 0 .5\n')
    N.select_policy(tmp_path)
    r=apply_thermal_policy(dict(charge=-1,mult=2,e_smd_eV=-10.,g_thermal_eV=2.,n_imag=0),tmp_path)
    assert r['legacy_g_thermal_eV']==2. and r['g_thermal_eV'] is None
    assert free_energy(r)==(None,'missing_validated_native_thermal')


def test_primary_native_requires_projected_partition(tmp_path):
    (tmp_path/'opt.xyz').write_text('3\n\nO 0 0 0\nH .7 0 .5\nH -.7 0 .5\n')
    p,r=record(tmp_path);p.write_text(json.dumps(r))
    assert N._load_native_surface_record(tmp_path,-1,2)
    assert N.load_native_thermal(tmp_path,-1,2) is None


def test_primary_native_uses_projected_values_and_retains_raw_evidence(tmp_path,monkeypatch):
    from redox.qm import projected_thermal as P
    from redox.qm.qrrho import PROTOCOL
    (tmp_path/'opt.xyz').write_text('3\n\nO 0 0 0\nH .7 0 .5\nH -.7 0 .5\n')
    p,r=record(tmp_path);p.write_text(json.dumps(r))
    post=tmp_path/'projected_thermal'/'post.json';post.parent.mkdir();post.write_text('{}')
    partition=dict(identity=dict(protocol=PROTOCOL),postprocessing_path=str(post),
                   nuclear_by_cutoff_eV={'25':1.6,'50':1.7,'100':1.8},frequencies_cm=[40.,90.,190.],
                   g_nuclear_qrrho_eV=1.7,g_spin_partition_eV=r['g_spin_partition_eV'],
                   g_thermal_total_eV=1.7+r['g_spin_partition_eV'])
    monkeypatch.setattr(P,'load_projected_thermal',lambda *args:partition)
    result=N.load_native_thermal(tmp_path,-1,2)
    assert result['g_nuclear_qrrho_eV']==1.7 and result['raw_native_g_thermal_total_eV']==r['g_thermal_total_eV']
    assert result['thermal_partition_protocol']==PROTOCOL
    assert result['thermal_partition_address']==p.stem+'-post'
    assert json.loads(p.read_text())==r


def test_certificate_allows_clean_checkout_and_rejects_changed_evidence(tmp_path,monkeypatch):
    import importlib,hashlib
    # The fixture stubs per-record control lookup; reload to test the real certificate admission.
    native=importlib.reload(N)
    certificate=json.loads((native.ROOT/'config/native_thermal_validation.json').read_text())
    (tmp_path/'config').mkdir();p=tmp_path/'config/native_thermal_validation.json';p.write_text(json.dumps(certificate))
    monkeypatch.setattr(native,'ROOT',tmp_path)
    monkeypatch.setattr(native,'installed_versions',lambda: {'pyscf':'2.14.0','gpu4pyscf':'1.8.1'})
    proof=native.implementation_controls();assert len(proof)==3
    assert [r['sha256'] for r in proof]==[r['source_sha256'] for r in certificate['controls']]
    source=tmp_path/certificate['controls'][0]['relative_path'];source.parent.mkdir(parents=True);source.write_text('{}')
    with pytest.raises(RuntimeError,match='differs'):native.implementation_controls()
    source.unlink();certificate['native_protocol']['grid_response']=False;p.write_text(json.dumps(certificate))
    with pytest.raises(RuntimeError,match='integrity'):native.implementation_controls()


def test_desired_spin_change_invalidates_existing_free_energy(tmp_path,monkeypatch):
    from redox.core import common
    from redox.qm import dft
    sd=tmp_path/'m/red2';sd.mkdir(parents=True)
    monkeypatch.setattr(common,'read_manifest',lambda:[dict(id='m',state='red2',mult_hint='3')])
    monkeypatch.setattr(dft,'_uma_mult',lambda *args:1)
    result=common.apply_thermal_policy(dict(charge=-2,mult=3,e_smd_eV=-10,g_thermal_eV=2),sd)
    assert result['desired_mult']==1
    assert common.free_energy(result)==(None,'wrong_spin')


def test_policy_selection_is_idempotent_and_updates_only_definition(tmp_path):
    N.select_policy(tmp_path);p=tmp_path/'thermal_policy.json'
    content=p.read_bytes();mtime=p.stat().st_mtime_ns
    N.select_policy(tmp_path)
    assert p.read_bytes()==content and p.stat().st_mtime_ns==mtime
    N.select_policy(tmp_path,'legacy_gfn2_gas_rrho')
    assert json.loads(p.read_text())['thermal_definition']=='legacy_gfn2_gas_rrho'
