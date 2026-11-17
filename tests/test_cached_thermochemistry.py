"""Portable arithmetic must reject tampering and identity/stoichiometry errors."""
import hashlib
import importlib.util
import json
from pathlib import Path
import shutil
import pytest

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('cached_reproduction',ROOT/'examples/cached_thermochemistry/reproduce.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)


@pytest.fixture
def cache(tmp_path):
    target=tmp_path/'moved';shutil.copytree(ROOT/'examples/cached_thermochemistry/data',target);return target


def test_relocated_no_quantum_dependencies(cache):
    result=module.reconstruct(cache)
    assert result['PMDI_first_reduction_V_vs_Fc']==pytest.approx(-1.0757908637970104,abs=1e-10)
    assert result['reactions'][1]['SI_deltaG_kJ_mol']==pytest.approx(-103.06882617,abs=1e-6)
    assert not result['quantum_calculations_run'] and not result['nuclear_partition_reintegrated']


def test_changed_raw_energy_refused(cache):
    manifest=json.loads((cache/'manifest.json').read_text());p=cache/manifest['states'][0]['files']['energy']
    p.write_text(p.read_text()+' ')
    with pytest.raises(ValueError,match='hash mismatch'):module.reconstruct(cache)


def test_rehashed_wrong_geometry_identity_refused(cache):
    manifest=json.loads((cache/'manifest.json').read_text());relative=manifest['states'][0]['files']['energy'];p=cache/relative
    raw=json.loads(p.read_text());raw['geom_sha1']='0'*16;p.write_text(json.dumps(raw))
    manifest['files_sha256'][relative]=hashlib.sha256(p.read_bytes()).hexdigest();(cache/'manifest.json').write_text(json.dumps(manifest))
    with pytest.raises(ValueError,match='does not match'):module.reconstruct(cache)


def test_unbalanced_reaction_refused(cache):
    p=cache/'manifest.json';manifest=json.loads(p.read_text());manifest['reactions'][0]['coefficients']['HCl']=2;p.write_text(json.dumps(manifest))
    with pytest.raises(ValueError,match='balance fails'):module.reconstruct(cache)
