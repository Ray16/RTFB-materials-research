import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace

import pytest


ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('cpu_ghost_pilot',ROOT/'scripts/pipeline/counterpoise_cpu_ghost_pilot.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)


def test_unfinished_own_term_cannot_release_ghost_pilot(tmp_path,monkeypatch):
    monkeypatch.setattr(module,'ROOT',tmp_path)
    monkeypatch.setattr(module.api,'load_module',lambda *args:SimpleNamespace(source_hashes=lambda *args:({},{})))
    own=tmp_path/'calcs'/'own';own.mkdir(parents=True);(own/'status.json').write_text(json.dumps({'status':'CPU_pilot_SCF_running_unqualified'}))
    with pytest.raises(ValueError,match='completed A-own'):
        module.inputs(tmp_path/'case',tmp_path/'geometry_review.json',own,tmp_path/'stationarity.json')


def test_prepared_ghost_geometry_has_full_basis_and_19_ghost_centers():
    core=module.api.load_module('test_large_ghost_protocol',ROOT/'scripts/validation/local_counterpoise.py')
    case=ROOT/'results/validation/local_counterpoise_2026-10-26_cpu_one_spacer/occ_2_4__pose1__cap109'
    geom=json.loads((case/'geometry.json').read_text());mols=core.molecules(geom['symbols'],geom['coordinates_A'],geom['A_indices'])
    ghost=mols['A_ghost'];assert ghost.nao_nr()==mols['AB'].nao_nr()==1026
    assert ghost.nelectron==84 and sum(ghost.atom_charges()==0)==19


def test_failed_own_stationarity_never_releases_ghost(tmp_path,monkeypatch):
    monkeypatch.setattr(module,'ROOT',tmp_path)
    monkeypatch.setattr(module.api,'load_module',lambda *args:SimpleNamespace(source_hashes=lambda *args:({},{})))
    control=tmp_path/'ghost_water.json';control.write_text(json.dumps({'status':'water_ghost_energy_and_orbital_stationarity_checked',
        'passed':True,'role':'A_ghost','orbital_gradient_norm_au':0.,'orbital_gradient_threshold_au':1e-7}))
    monkeypatch.setattr(module,'GHOST_CONTROL',control)
    own=tmp_path/'calcs'/'own';(own/'A_own').mkdir(parents=True)
    state={'status':'CPU_pilot_converged_independent_review_pending','converged':True,'qualified':False,'source_sha256':{},'E_Ha':-1.}
    (own/'status.json').write_text(json.dumps(state));(own/'A_own/energy.json').write_text(json.dumps(state));(own/'launch.json').write_text(json.dumps({'source_sha256':{}}))
    gradient=tmp_path/'gradient.json';gradient.write_text(json.dumps({'status':'CPU_pilot_energy_and_orbital_stationarity_checked',
        'passed':False,'full_counterpoise_computed':False,'solvent_counterpoise_computed':False,'reference_energy_Ha':-1.}))
    with pytest.raises(ValueError,match='stationary A-own'):
        module.inputs(tmp_path/'case',tmp_path/'geometry.json',own,gradient)


def test_stationarity_must_bind_the_exact_checkpoint(tmp_path,monkeypatch):
    monkeypatch.setattr(module,'ROOT',tmp_path)
    monkeypatch.setattr(module.api,'load_module',lambda *args:SimpleNamespace(source_hashes=lambda *args:({},{})))
    monkeypatch.setattr(module.api,'verify_sources',lambda *args:None)
    control=tmp_path/'ghost_water.json';control.write_text(json.dumps({'status':'water_ghost_energy_and_orbital_stationarity_checked',
        'passed':True,'role':'A_ghost','orbital_gradient_norm_au':0.,'orbital_gradient_threshold_au':1e-7,'source_sha256':{}}))
    monkeypatch.setattr(module,'GHOST_CONTROL',control)
    own=tmp_path/'calcs'/'own';(own/'A_own').mkdir(parents=True)
    state={'status':'CPU_pilot_converged_independent_review_pending','converged':True,'qualified':False,'source_sha256':{},'E_Ha':-1.}
    (own/'status.json').write_text(json.dumps(state));(own/'A_own/energy.json').write_text(json.dumps(state));(own/'launch.json').write_text(json.dumps({'source_sha256':{}}))
    for path in (own/'A_own/SCF.h5',own/'grid_coordinates_AU.npy',own/'grid_weights.npy'):path.write_bytes(b'fixture')
    gradient=tmp_path/'gradient.json';gradient.write_text(json.dumps({'status':'CPU_pilot_energy_and_orbital_stationarity_checked',
        'passed':True,'full_counterpoise_computed':False,'solvent_counterpoise_computed':False,'reference_energy_Ha':-1.,'source_sha256':{}}))
    with pytest.raises(ValueError,match='exact A-own checkpoint'):
        module.inputs(tmp_path/'case',tmp_path/'geometry.json',own,gradient)


def test_ghost_and_own_numerical_environments_must_match(tmp_path,monkeypatch):
    monkeypatch.setattr(module,'ROOT',tmp_path)
    monkeypatch.setattr(module,'inputs',lambda *args:({},{},{}))
    own=tmp_path/'calcs'/'own';own.mkdir(parents=True)
    (own/'launch.json').write_text(json.dumps({'versions':{'pyscf':'wrong','numpy':'wrong'}}))
    with pytest.raises(ValueError,match='environments differ'):
        module.launch(tmp_path/'case',tmp_path/'geometry.json',own,tmp_path/'gradient.json',tmp_path/'calcs'/'ghost',8)
