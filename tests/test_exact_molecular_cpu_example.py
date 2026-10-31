import importlib.util
import json
from pathlib import Path
import shutil

import numpy as np
import pytest

ROOT=Path(__file__).parents[1]
spec=importlib.util.spec_from_file_location('portable_molecular_replay',ROOT/'examples/exact_viologen_replay/replay.py')
reader=importlib.util.module_from_spec(spec);spec.loader.exec_module(reader)


@pytest.fixture
def package(tmp_path):
    target=tmp_path/'package';shutil.copytree(ROOT/'examples/exact_viologen_replay/data',target);return target


def test_ten_cached_states_reconstruct_exact_redox_and_grafting_without_writing_inputs(package):
    before={str(p.relative_to(package)):reader.sha(p) for p in package.rglob('*') if p.is_file()}
    result=reader.replay(package)
    assert result['states']==10 and result['Hessians_reconstructed']==9 and result['FD_states_reconstructed']==6
    assert result['by_cutoff']['50']['E1_V_vs_Fc']==pytest.approx(-.6077915087)
    assert result['by_cutoff']['50']['graft_precursor_2_kJ_mol']==pytest.approx(-11.49332965)
    assert before=={str(p.relative_to(package)):reader.sha(p) for p in package.rglob('*') if p.is_file()}
    assert not result['material_admission'] and not result['quantum_calculations_run']


def test_numpy_scalars_are_valid_numbers_and_boolean_metadata_cannot_replace_them():
    reader.close(np.float64(1),1.)
    with pytest.raises(ValueError):reader.close(True,1.)


@pytest.mark.parametrize('change',['hash','charge','scope','identity','protocol','nan_frequency','nuclear','FD_gradient','FD_settings','stoichiometry','selection','unknown_kind','reported_bad_FD'])
def test_changed_physics_records_or_scope_cannot_pass_after_payload_rehash(package,change):
    p=package/'manifest.json';manifest=json.loads(p.read_text());state=manifest['states'][0];modified=None
    if change=='hash':(package/state['files']['geometry']).write_text('changed')
    if change=='charge':state['charge']=1
    if change=='scope':manifest['full_historical_source_closure_revalidated']=True
    if change=='identity':manifest['product_inchikey']='different'
    if change=='protocol':manifest['active_protocol']['basis']='def2-tzvp'
    if change in ('nan_frequency','nuclear'):
        modified=state['files']['partition'];path=package/modified;d=json.loads(path.read_text())
        if change=='nan_frequency':d['frequencies_cm'][0]=float('nan')
        else:d['nuclear_by_cutoff_eV']['50']+=.01
        path.write_text(json.dumps(d))
    if change=='FD_gradient':
        modified=state['files']['fd_arrays'];path=package/modified;data=np.load(path,allow_pickle=False);coords=data['coordinates_Bohr'];g=data['gradient_Ha_Bohr'].copy();g[1,0,0]+=.01;np.savez_compressed(path,coordinates_Bohr=coords,gradient_Ha_Bohr=g)
    if change=='FD_settings':
        modified=state['files']['fd_records'];path=package/modified;d=json.loads(path.read_text());d['points'][0]['settings']['charge']=0;path.write_text(json.dumps(d))
    if change=='stoichiometry':manifest['reactions'][2]['coefficients'].pop('chloride')
    if change=='selection':manifest['selected_precursor_id']='precursor_0'
    if change=='unknown_kind':manifest['reactions'][0]['kind']='other'
    if change=='reported_bad_FD':
        modified=state['files']['native'];path=package/modified;d=json.loads(path.read_text());d['derivative_checks'][0]['relative_error']=.02;path.write_text(json.dumps(d))
    if modified:manifest['files_sha256'][modified]=reader.sha(package/modified)
    p.write_text(json.dumps(manifest))
    with pytest.raises(ValueError):reader.replay(package)


@pytest.mark.parametrize('relative',['/absolute','../escape','a/../../escape'])
def test_package_paths_cannot_escape_the_supplied_root(tmp_path,relative):
    with pytest.raises(ValueError):reader.local(tmp_path,relative)
