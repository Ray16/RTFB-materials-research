"""Equilibrium spin populations need bound component evidence and no repeated degeneracy."""
import copy
import json
import math

import pytest

from redox.core.common import free_energy
from redox.validation import spin_ensemble as ensemble


def test_equal_electronic_and_nuclear_terms_include_triplet_degeneracy_once():
    kt=8.617333262145e-5*298.15
    result=ensemble.equilibrium(-10.,-10.-kt*math.log(3))
    assert result['alternative_population']==pytest.approx(.75)
    assert result['free_energy_eV']==pytest.approx(-10.-kt*math.log(4))


@pytest.mark.parametrize('gap',[-1000.,1000.])
def test_extreme_gaps_are_stable_and_choose_the_lower_surface(gap):
    result=ensemble.equilibrium(-1e6,-1e6+gap)
    assert math.isfinite(result['free_energy_eV'])
    assert result['alternative_population']==float(gap<0)
    assert result['free_energy_eV']==min(-1e6,-1e6+gap)


@pytest.mark.parametrize('value',[float('nan'),float('inf'),True])
def test_nonphysical_component_input_is_refused(value):
    with pytest.raises(ValueError):ensemble.equilibrium(value,0.)


@pytest.fixture
def selection(tmp_path,monkeypatch):
    artifact=tmp_path/'component.json';artifact.write_text('{}')
    components=[dict(role='chosen',free_energy_eV=-10.,electronic_eV=-11.,thermal_eV=1.),
                dict(role='alternative',free_energy_eV=-9.98,electronic_eV=-10.9,thermal_eV=.92)]
    current=dict(chosen_identity={'charge':-2,'mult':1},alternative_seed_identity={'geometry':'old seed'},
                 helper='/helper',components=components,
                 native_partition_protocol={'temperature_K':298.15},
                 evidence_sha256={str(artifact):ensemble.digest(artifact)})
    monkeypatch.setattr(ensemble,'_current_components',lambda *args:copy.deepcopy(current))
    policy=ensemble.prepare_policy(tmp_path,'/helper')
    record=dict(charge=-2,mult=1,e_smd_eV=-11.,g_thermal_eV=1.,n_imag=0,thermal_qc='ok',
                thermal_definition='native_solution_dft_qrrho_spin',
                native_thermal_record={'g_thermal_total_eV':1.})
    return tmp_path,artifact,current,policy,record


def test_selected_mixture_changes_primary_G_and_preserves_single_surface_for_gaps(selection):
    path,_,_,policy,record=selection
    (path/ensemble.POLICY_NAME).write_text(json.dumps(policy))
    output=ensemble.apply_selected_ensemble(record,path)
    assert output['spin_ensemble_status']=='applied'
    assert free_energy(output)[0]==pytest.approx(policy['equilibrium']['free_energy_eV'])
    assert output['spin_ensemble_selected_surface_G_eV']==-10.
    assert output['single_surface_kinetics_validated'] is False
    assert record['g_thermal_eV']==1.  # no mutation of the selected surface


@pytest.mark.parametrize('change',['component','saved_population','evidence','temperature','kinetic_claim'])
def test_stale_or_tampered_policy_withholds_G(selection,change):
    path,artifact,current,policy,record=selection
    if change=='component':current['components'][1]['free_energy_eV']+=.01
    elif change=='saved_population':policy['equilibrium']['alternative_population']+=.01
    elif change=='evidence':artifact.write_text('{"changed":true}')
    elif change=='temperature':policy['temperature_K']=300.
    else:policy['single_surface_kinetics_validated']=True
    (path/ensemble.POLICY_NAME).write_text(json.dumps(policy))
    output=ensemble.apply_selected_ensemble(record,path)
    assert output['spin_ensemble_status']=='invalid'
    assert free_energy(output)[0] is None


def test_absent_policy_cannot_inherit_forged_ensemble_status(selection):
    path,_,_,_,record=selection
    record.update(spin_ensemble_status='applied',spin_ensemble_free_energy_eV=-100.)
    output=ensemble.apply_selected_ensemble(record,path)
    assert output['spin_ensemble_status']=='not_selected'
    assert 'spin_ensemble_free_energy_eV' not in output
    assert free_energy(output)[0]==-10.


def test_selected_policy_never_recovers_an_unqualified_native_component(selection):
    path,_,_,policy,record=selection
    (path/ensemble.POLICY_NAME).write_text(json.dumps(policy))
    record['native_thermal_record']=None
    output=ensemble.apply_selected_ensemble(record,path)
    assert free_energy(output)[0] is None


def test_ensemble_admission_requires_the_applied_primary_value(selection):
    path,_,_,policy,record=selection
    (path/ensemble.POLICY_NAME).write_text(json.dumps(policy))
    output=ensemble.apply_selected_ensemble(record,path)
    assert ensemble.admitted_ensemble(output,path) is not None
    output['g_thermal_eV']=record['g_thermal_eV']
    assert ensemble.admitted_ensemble(output,path) is None
    output['g_thermal_eV']=float('nan')
    assert ensemble.admitted_ensemble(output,path) is None


@pytest.fixture
def independent_case(tmp_path,monkeypatch):
    import hashlib
    import importlib.util
    from redox.core import common,protocol
    from pathlib import Path
    root=tmp_path;chosen=root/'calcs/dft/m/neu';helper=root/'calcs/spincheck/m/neu_m3'
    sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    native_protocol={'schema':1,'model':'fixture native contract'}
    partition_protocol={'schema':1,'temperature_K':298.15}
    evidence={};components=[];thermal={};partitions={}
    identity=None;seed_identity=None
    for role,folder,m,e,g in [('chosen',chosen,1,-11.,1.),('alternative',helper,3,-10.9,.92)]:
        folder.mkdir(parents=True)
        xyz=folder/'opt.xyz';xyz.write_text('2\n\nH 0 0 0\nH 0 0 .74\n')
        raw=dict(charge=0,mult=m,optimized=True,opt_xc='r2scan',opt_basis='def2-svp',
                 opt_disp='d4',solvent='acetonitrile')
        sp=protocol.record_path(folder,0,m);sp.parent.mkdir()
        sp.write_text(json.dumps(dict(protocol=protocol.ACTIVE_SP,protocol_hash=sp.stem,
            geom_sha1=protocol.geom_sha1(xyz),charge=0,mult=m,converged_gas=True,converged_smd=True,
            e_smd_eV=e,e_gas_eV=e,smd_spin_contam=0.)))
        if role=='chosen':
            identity=dict(id='m',state='neu',charge=0,mult=1,geometry_sha256=sha(xyz),
                smiles_sha256=hashlib.sha256(b'[H][H]').hexdigest(),sp_record_sha256=sha(sp),protocol=protocol.ACTIVE_SP)
        else:
            seed_identity=dict(identity,geometry_sha256='earlier geometry',sp_record_sha256='earlier energy')
            raw['spin_source_identity']=seed_identity
        rj=folder/'result.json';rj.write_text(json.dumps(raw))
        tp=folder/'thermal_policy.json';tp.write_text(json.dumps(dict(
            thermal_definition='native_solution_dft_qrrho_spin',native_protocol=native_protocol)))
        ident=dict(geometry_sha256=sha(xyz),charge=0,mult=m,protocol=native_protocol)
        address=hashlib.sha256(json.dumps(ident,sort_keys=True).encode()).hexdigest()[:16]
        native=folder/'native_thermal'/f'{address}.json';native.parent.mkdir();native.write_text('{}')
        h=native.with_suffix('')/'hessian_cartesian.npy';h.parent.mkdir();h.write_bytes(b'fixture')
        post=folder/'projected_thermal/p.json';post.parent.mkdir();artifact=post.with_suffix('')/'frequencies.json'
        artifact.parent.mkdir();artifact.write_text('{}')
        post.write_text(json.dumps({'artifacts_sha256':{str(artifact):sha(artifact)}}))
        paths=[rj,xyz,tp,sp,native,h,post,artifact];evidence.update({str(p):sha(p) for p in paths})
        components.append(dict(role=role,directory=str(folder),charge=0,mult=m,electronic_eV=e,
            thermal_eV=g,free_energy_eV=e+g,geometry_sha256=sha(xyz),energy_path=str(sp),
            native_path=str(native),partition_path=str(post)))
        thermal[str(folder)]=(g,0,'native_solution_dft_qrrho_spin')
        partitions[str(folder)]=dict(protocol=partition_protocol,path=str(post))
    current=dict(chosen_identity=identity,alternative_seed_identity=seed_identity,helper=str(helper),
        components=components,native_partition_protocol=partition_protocol,evidence_sha256=evidence)
    monkeypatch.setattr(ensemble,'_current_components',lambda *args:copy.deepcopy(current))
    monkeypatch.setattr(common,'smiles_to_xyz_map',lambda *args:[0,1])
    policy=ensemble.prepare_policy(chosen,helper)
    spec=importlib.util.spec_from_file_location('independent_ensemble_test',common.ROOT/'scripts/validation/audit/spin_ensemble.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    kwargs=dict(root=root,manifest=[dict(id='m',state='neu',charge=0,smiles='[H][H]')],
        native_protocol=native_protocol,surface_reader=lambda p,r:thermal[str(p)],
        partition_reader=lambda p,r:partitions[str(p)],spin_decision=lambda *args:True)
    return module,chosen,helper,policy,kwargs,thermal


def test_independent_reader_retains_original_seed_and_replays_both_own_minima(independent_case):
    module,chosen,_,policy,kwargs,_=independent_case
    result=module.read_ensemble(chosen,policy,**kwargs,strict=True)
    assert result['free_energy_eV']==pytest.approx(policy['equilibrium']['free_energy_eV'])
    assert policy['alternative_seed_identity']['geometry_sha256']=='earlier geometry'


@pytest.mark.parametrize('change',['component','seed_chemistry','native_missing','policy_artifact','wrong_spin'])
def test_independent_reader_refuses_unqualified_or_changed_components(independent_case,change):
    module,chosen,helper,policy,kwargs,thermal=independent_case
    if change=='component':policy['components'][1]['free_energy_eV']+=.001
    elif change=='seed_chemistry':
        p=helper/'result.json';raw=json.loads(p.read_text());raw['spin_source_identity']['smiles_sha256']='different'
        p.write_text(json.dumps(raw))
    elif change=='native_missing':thermal[str(helper)]=(None,None,'native_solution_dft_qrrho_spin')
    elif change=='policy_artifact':(helper/'projected_thermal/p/frequencies.json').write_text('{"changed":true}')
    else:kwargs['spin_decision']=lambda *args:False
    assert module.read_ensemble(chosen,policy,**kwargs) is None


def test_equilibrium_population_keeps_gas_lambda_as_a_flagged_diagnostic(tmp_path,monkeypatch):
    from redox.properties import reorg
    ox=tmp_path/'ox.json';red=tmp_path/'red.json'
    ox.write_text(json.dumps({'e_gas_eV':-9.8}))
    red.write_text(json.dumps({'e_gas_eV':-10.7,'gas_homo_eV':-1.}))
    states={'neu':{'e_gas_eV':-10.},'red1':{'e_gas_eV':-11.,'gas_homo_eV':-1.,
                                         'spin_ensemble_status':'applied'}}
    monkeypatch.setattr(reorg,'_res',lambda gid,st:states[st])
    monkeypatch.setattr(reorg,'_cache_path',lambda gid,st,*args:ox if st=='neu' else red)
    monkeypatch.setattr(reorg,'_cache_ok',lambda *args:True)
    monkeypatch.setattr(reorg,'_geom_hash',lambda *args:'geometry')
    monkeypatch.setattr(reorg,'_heavy_rmsd',lambda *args:.1)
    monkeypatch.setattr(reorg,'_sym_heavy_rmsd',lambda *args,**kw:.1)
    row=reorg.lambda_for_couple('m',('neu',0,1),('red1',-1,2))
    assert row['lambda_i_eV']==.5
    assert row['flag']=='spin_ensemble_kinetics_unvalidated'


def test_selection_requires_the_shared_reader_hook_before_any_write(selection,monkeypatch):
    from redox.core import common
    path,_,_,policy,_=selection
    monkeypatch.delattr(common,'_apply_single_surface_thermal_policy',raising=False)
    with pytest.raises(RuntimeError,match='not activated'):ensemble.select_policy(path,policy)
    assert not (path/ensemble.POLICY_NAME).exists()


@pytest.fixture
def selectable(selection,monkeypatch):
    from redox.core import common
    path,_,_,policy,record=selection
    monkeypatch.setattr(common,'ROOT',path)
    monkeypatch.setattr(common,'_apply_single_surface_thermal_policy',lambda *args:None,raising=False)
    monkeypatch.setattr(common,'read_result',lambda *args:ensemble.apply_selected_ensemble(record,path))
    monkeypatch.setattr(ensemble,'_independent_thermal_inputs',lambda *args:
                        (policy['equilibrium']['free_energy_eV']-record['e_smd_eV'],0,'native_solution_dft_qrrho_spin'))
    return selection


def test_selection_commits_only_after_both_readers_admit(selectable):
    path,_,_,policy,_=selectable
    result=ensemble.select_policy(path,policy)
    assert result['status']=='equilibrium_policy_selected_and_both_readers_verified'
    assert (path/ensemble.POLICY_NAME).exists()


@pytest.mark.parametrize('previous',[None,b'{"old":"selection"}',b'\xffprior bytes'])
def test_failed_independent_selection_restores_exact_previous_bytes(selectable,monkeypatch,previous):
    path,_,_,policy,_=selectable
    target=path/ensemble.POLICY_NAME
    if previous is not None:target.write_bytes(previous)
    monkeypatch.setattr(ensemble,'_independent_thermal_inputs',lambda *args:(None,None,'native_solution_dft_qrrho_spin'))
    with pytest.raises(RuntimeError,match='Independent reader'):ensemble.select_policy(path,policy)
    assert target.read_bytes()==previous if previous is not None else not target.exists()
    assert list((path/'spin_ensemble_history').glob('*.failed.json'))
