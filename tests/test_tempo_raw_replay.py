import importlib.util
import json
import os
from pathlib import Path
import shutil
import hashlib

import pytest


@pytest.fixture
def package(tmp_path):
    source=os.environ.get('REDOX_TEMPO_REPLAY_PACKAGE')
    if not source:pytest.skip('Set REDOX_TEMPO_REPLAY_PACKAGE to the extracted local addendum')
    source=Path(source)
    if not (source/'manifest.json').is_file():pytest.skip('TEMPO raw package is unavailable')
    target=tmp_path/'package';shutil.copytree(source,target);return target


def load(package):
    spec=importlib.util.spec_from_file_location('tempo_replay_test',package/'reproduce.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module


def rewrite_manifest(package,manifest):
    (package/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')


def rehash(package,manifest,relative):
    manifest['files_sha256'][relative]=hashlib.sha256((package/relative).read_bytes()).hexdigest()


def test_reconstructs_four_states_and_three_cutoffs(package):
    result=load(package).replay(package,full_native=True)
    assert result['oxidation_V_vs_Fc_by_cutoff']['50']==pytest.approx(.5787750407616841,abs=1e-9)
    assert result['native_reconstruction']['Hessians_reconstructed']==4
    assert not result['experimental_accuracy_inferred'] and not result['material_eligibility_inferred']


@pytest.mark.parametrize('mutation',['scope','duplicate','path_escape','energy','hessian'])
def test_rejects_changed_or_promoted_evidence(package,mutation):
    manifest=json.loads((package/'manifest.json').read_text())
    if mutation=='scope':manifest['material_eligibility_inferred']=True
    if mutation=='duplicate':manifest['states'][1]['id']=manifest['states'][0]['id']
    if mutation=='path_escape':manifest['states'][0]['files']['geometry']='../outside.xyz'
    if mutation=='energy':
        relative=manifest['states'][0]['files']['energy'];p=package/relative;r=json.loads(p.read_text());r['e_smd_eV']+=.1;p.write_text(json.dumps(r));rehash(package,manifest,relative)
    if mutation=='hessian':
        import numpy as np
        relative=manifest['states'][0]['files']['hessian'];p=package/relative;arr=np.load(p,allow_pickle=False);arr[0,0]+=.01;np.save(p,arr);rehash(package,manifest,relative)
    rewrite_manifest(package,manifest)
    with pytest.raises((ValueError,FileNotFoundError)):
        load(package).replay(package,full_native=(mutation=='hessian'))
