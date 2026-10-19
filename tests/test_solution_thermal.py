import copy
import json
from redox.qm.thermal import thermal_identity,thermal_hash,spectrum_values,stationary_verdict,paired_potential,load_solution_thermal_record,spin_degeneracy_eV


def stable():
    return dict(returncode=0,optimization_converged=True,connectivity_intact=True,
        frequencies_cm=[0,0,0,0,0,0,30,400],n_imag=0,g_thermal_eV=2.5)


def test_rrho_count_cannot_hide_negative_mode():
    rec=stable();rec['frequencies_cm'].append(-11.62)
    ok,issues=stationary_verdict(rec)
    assert not ok and 'unstable_mode' in issues


def test_thermal_source_and_spin_define_identity(tmp_path):
    p=tmp_path/'geom.xyz';p.write_text('1\n\nH 0 0 0\n')
    singlet=thermal_identity(p,-2,1);triplet=thermal_identity(p,-2,3)
    assert thermal_hash(singlet)!=thermal_hash(triplet)
    rec=stable();rec['identity']=triplet
    assert 'source_or_protocol_mismatch' in stationary_verdict(rec,singlet)[1]
    p.write_text('1\n\nH 0.1 0 0\n')
    assert thermal_hash(thermal_identity(p,-2,1))!=thermal_hash(singlet)


def test_decomposition_and_failed_execution_block_model_minimum():
    for key,value,issue in [('connectivity_intact',False,'changed_chemistry'),('returncode',1,'failed_execution'),('optimization_converged',False,'nonstationary')]:
        rec=stable();rec[key]=value
        assert issue in stationary_verdict(rec)[1]


def test_matched_fc_cancels_uniform_standard_state_shift():
    original=paired_potential(8,3,1,2,7,3,4,5)
    shifted=paired_potential(8,3,1.083,2.083,7,3,4.083,5.083)
    assert abs(original-shifted)<1e-12
    assert paired_potential(8,3,1,None,7,3,4,5) is None


def test_spectrum_parser_keeps_reported_imaginary_value():
    assert spectrum_values('$vibrational spectrum\n1 -0.00 0 -\n7 a -11.62 2 YES\n')==[-0.,-11.62]


def test_stationarity_does_not_imply_validated_transfer(tmp_path):
    geom=tmp_path/'opt.xyz';geom.write_text('2\n\nH 0 0 0\nH 0 0 0.7\n')
    ident=thermal_identity(geom,0,1);rec=stable();rec['identity']=ident
    rec['frequencies_cm']=[0,0,0,0,0,400]
    p=tmp_path/'selected.json';p.write_text(json.dumps(dict(identity=ident,selected=rec,g_thermal_by_cutoff_eV={'50':2.5})))
    assert load_solution_thermal_record(p,geom,0,1)
    assert load_solution_thermal_record(p,geom,0,1,require_transfer_validation=True) is None
    assert load_solution_thermal_record(p,geom,0,3) is None
    geom.write_text('2\n\nH 0.1 0 0\nH 0 0 0.7\n')
    assert load_solution_thermal_record(p,geom,0,1) is None


def test_spin_partition_reference_cancels_second_quinone_wave_only():
    singlet=spin_degeneracy_eV(1);doublet=spin_degeneracy_eV(2)
    assert singlet==0 and abs(doublet+0.017808737388)<1e-7
    first=(singlet-doublet)-(doublet-singlet)
    second=(doublet-singlet)-(doublet-singlet)
    assert abs(first-0.035617474776)<1e-7 and second==0
