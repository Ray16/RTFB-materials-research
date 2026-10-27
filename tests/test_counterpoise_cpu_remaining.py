import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace

import pytest


ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('cpu_remaining',ROOT/'scripts/pipeline/counterpoise_cpu_remaining.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)


def test_prepared_large_geometry_matches_all_remaining_role_identities():
    core=module.api.load_module('remaining_geometry_protocol_test',ROOT/'scripts/validation/local_counterpoise.py')
    case=ROOT/'results/validation/local_counterpoise_2026-10-26_cpu_one_spacer/occ_2_4__pose1__cap109'
    geometry=json.loads((case/'geometry.json').read_text());mols=core.molecules(geometry['symbols'],geometry['coordinates_A'],geometry['A_indices'])
    for role,(electrons,basis,ghosts) in module.EXPECTED.items():
        mol=mols[role]
        assert mol.nelectron==electrons and mol.nao_nr()==basis and sum(mol.atom_charges()==0)==ghosts


def test_incomplete_reference_cannot_release_remaining_terms(tmp_path,monkeypatch):
    monkeypatch.setattr(module,'ROOT',tmp_path)
    monkeypatch.setattr(module.api,'load_module',lambda *args:SimpleNamespace(source_hashes=lambda *args:({},{})))
    monkeypatch.setattr(module.api,'verify_sources',lambda *args:None)
    control=tmp_path/'control.json';control.write_text(json.dumps({'status':'remaining_water_counterpoise_terms_stationary',
        'passed':True,'source_sha256':{},'terms':[{'role':r,'passed':True} for r in module.EXPECTED]}));monkeypatch.setattr(module,'CONTROL',control)
    own=tmp_path/'own';own.mkdir();(own/'status.json').write_text(json.dumps({'status':'CPU_pilot_SCF_running_unqualified'}))
    with pytest.raises(ValueError,match='Both matched A terms'):
        module.inputs(tmp_path/'case',tmp_path/'review',own,tmp_path/'own_gradient',tmp_path/'ghost',tmp_path/'ghost_gradient')


def test_unknown_role_is_rejected_before_launch(tmp_path,monkeypatch):
    monkeypatch.setattr(module,'ROOT',tmp_path)
    with pytest.raises(ValueError,match='Supported term'):
        module.launch('A_own',tmp_path/'case',tmp_path/'review',tmp_path/'own',tmp_path/'own_gradient',
                      tmp_path/'ghost',tmp_path/'ghost_gradient',tmp_path/'calcs'/'task',8)


def test_failed_reference_stationarity_cannot_release_remaining_terms(tmp_path,monkeypatch):
    monkeypatch.setattr(module,'ROOT',tmp_path)
    monkeypatch.setattr(module.api,'load_module',lambda *args:SimpleNamespace(source_hashes=lambda *args:({},{})))
    monkeypatch.setattr(module.api,'verify_sources',lambda *args:None)
    control=tmp_path/'control.json';control.write_text(json.dumps({'status':'remaining_water_counterpoise_terms_stationary',
        'passed':True,'source_sha256':{},'terms':[{'role':r,'passed':True} for r in module.EXPECTED]}));monkeypatch.setattr(module,'CONTROL',control)
    own=tmp_path/'own';(own/'A_own').mkdir(parents=True)
    state={'status':'CPU_pilot_converged_independent_review_pending','converged':True,'qualified':False,'source_sha256':{},'E_Ha':-1.}
    (own/'status.json').write_text(json.dumps(state));(own/'A_own/energy.json').write_text(json.dumps(state));(own/'launch.json').write_text(json.dumps({'source_sha256':{}}))
    gradient=tmp_path/'gradient.json';gradient.write_text(json.dumps({'status':'CPU_pilot_energy_and_orbital_stationarity_checked',
        'passed':False,'full_counterpoise_computed':False,'solvent_counterpoise_computed':False,'reference_energy_Ha':-1.}))
    with pytest.raises(ValueError,match='stationarity gate'):
        module.inputs(tmp_path/'case',tmp_path/'review',own,gradient,tmp_path/'ghost',tmp_path/'ghost_gradient')
