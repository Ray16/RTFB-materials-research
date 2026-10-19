import math,json
import pytest
from redox.validation import spins


def test_ensemble_counts_singlet_triplet_spin_weights_once():
    kt=8.617333262145e-5*298.15
    result=spins.two_state_ensemble(0.,-kt*math.log(3))
    assert result['ensemble_free_energy_eV']==pytest.approx(-kt*math.log(4))
    assert result['alternative_population']==pytest.approx(.75)


def test_ensemble_stable_for_inaccessible_alternative_and_bad_temperature():
    result=spins.two_state_ensemble(-100000.,-99000.)
    assert result['alternative_population']==0 and result['ensemble_correction_from_electronic_selected_eV']==0
    with pytest.raises(ValueError):spins.two_state_ensemble(1,2,0)


def test_source_addressed_spin_report_rejects_changed_chosen_state(tmp_path,monkeypatch):
    monkeypatch.setattr(spins,'DFT',tmp_path)
    ident=dict(id='m',state='red2',geometry_sha256='old',charge=-2,mult=1)
    monkeypatch.setattr(spins,'snapshot',lambda *args:ident)
    path=spins.comparison_path('m','red2',ident);path.parent.mkdir(parents=True)
    path.write_text(json.dumps(dict(chosen_identity=ident,status='spin_gap_verified',electronic_gap_eV=.8)))
    assert spins.read_comparison('m','red2')
    ident['geometry_sha256']='new'
    assert spins.read_comparison('m','red2') is None


@pytest.mark.parametrize('chosen,alternative',[(1,3),(3,1),(2,4),(4,2)])
def test_compared_spins_preserve_electron_parity(chosen,alternative):
    assert spins.alternate_mult(chosen)==alternative


@pytest.fixture
def reproduced_gap():
    old=dict(id='m',state='red2',charge=-2,mult=1,geometry_sha256='g',protocol={'xc':'same'},sp_record_sha256='old')
    new=dict(old,sp_record_sha256='new')
    energy=dict(e_smd_eV=-100.,e_gas_eV=-99.)
    chosen=dict(energy,e_smd_eV=-99.9999)
    alt=dict(e_smd_eV=-99.3,smd_spin_contam=.01,mult=3,charge=-2)
    report=dict(status='spin_gap_verified',chosen_identity=old,chosen_mult=1,alternative_mult=3,
                electronic_gap_eV=.7,thermochemical_gap_validated=False)
    return [old,new,energy,chosen,alt,report]


def test_reproduced_sp_refresh_recalculates_gap_and_preserves_original_report(reproduced_gap):
    args=reproduced_gap;before=json.dumps(args[-1],sort_keys=True)
    assert spins.refreshed_electronic_gap(*args)==pytest.approx(.6999)
    assert json.dumps(args[-1],sort_keys=True)==before


@pytest.mark.parametrize('index,change',[(1,{'geometry_sha256':'changed'}),(1,{'protocol':{'xc':'different'}}),
    (3,{'e_smd_eV':-99.98}),(4,{'smd_spin_contam':.3}),(4,{'e_smd_eV':-99.2}),
    (5,{'status':'running'}),(5,{'thermochemical_gap_validated':True}),
    (4,{'mult':1})])
def test_refresh_refuses_changed_science_or_incomplete_evidence(reproduced_gap,index,change):
    reproduced_gap[index].update(change)
    with pytest.raises(ValueError):spins.refreshed_electronic_gap(*reproduced_gap)


def test_refresh_does_not_admit_a_near_degenerate_comparison(reproduced_gap):
    args=reproduced_gap;args[4]['e_smd_eV']=-99.9;args[5]['electronic_gap_eV']=.1
    with pytest.raises(ValueError,match='thermochemistry'):spins.refreshed_electronic_gap(*args)
