import hashlib
import json
from pathlib import Path
import shutil
import types

import numpy as np
import pytest

ROOT=Path(__file__).resolve().parents[1]
PACKAGE=ROOT/'examples/exact_viologen_corrected_replay/data'


def replay(path):
    source=path/'replay.py';module=types.ModuleType('checked_corrected_primary_example');module.__file__=str(source)
    exec(compile(source.read_bytes(),str(source),'exec'),module.__dict__)
    return module.replay(path)


def test_all_nine_new_native_gradient_panels_replay_without_project_runtime():
    r=replay(PACKAGE)
    assert (r['states'],r['Hessians_reconstructed'],r['FD_states_reconstructed'])==(10,9,9)
    assert r['selected_precursor_id']=='precursor_2' and not r['quantum_calculations_run']
    assert not r['material_admission'] and not r['full_historical_source_closure_revalidated']


def test_rehashed_changed_Fc_gradient_still_fails_physical_reconstruction(tmp_path):
    package=tmp_path/'data';shutil.copytree(PACKAGE,package)
    p=package/'states/fc_ox/fd_arrays.npz'
    with np.load(p) as a:coords=a['coordinates_Bohr'];gradients=a['gradient_Ha_Bohr']
    gradients[1,0,0]+=.001
    np.savez_compressed(p,coordinates_Bohr=coords,gradient_Ha_Bohr=gradients)
    m=package/'manifest.json';data=json.loads(m.read_text());data['files_sha256'][str(p.relative_to(package))]=hashlib.sha256(p.read_bytes()).hexdigest();m.write_text(json.dumps(data))
    with pytest.raises(ValueError,match='Cached quantity|Saved FD'):replay(package)


def test_numerical_cache_cannot_promote_material_or_sampling_scope(tmp_path):
    package=tmp_path/'data';shutil.copytree(PACKAGE,package);p=package/'manifest.json';m=json.loads(p.read_text());m['material_admission']=True;p.write_text(json.dumps(m))
    with pytest.raises(ValueError,match='scope'):replay(package)
