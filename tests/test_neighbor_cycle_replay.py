"""Integration checks against the independently packaged fourteen raw states.

Set REDOX_NEIGHBOR_REPLAY_PACKAGE to an extracted archive; these tests do not
read original quantum directories or require quantum/scientific packages.
"""
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import shutil

import pytest

ROOT=Path(__file__).parents[1]
spec=importlib.util.spec_from_file_location('neighbor_cached_replay_test',ROOT/'examples/neighbor_cycle_replay/reproduce.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)


@pytest.fixture
def package(tmp_path):
    source=os.environ.get('REDOX_NEIGHBOR_REPLAY_PACKAGE')
    if not source:pytest.skip('Extracted raw-state package must be supplied explicitly')
    target=tmp_path/'package';shutil.copytree(source,target);return target


def write(path,value):path.write_text(json.dumps(value,indent=2)+'\n')


def rehash(root,relative,mutate):
    path=root/relative;value=json.loads(path.read_text());mutate(value);write(path,value)
    manifest=json.loads((root/'manifest.json').read_text())
    manifest['files_sha256'][relative]=hashlib.sha256(path.read_bytes()).hexdigest();write(root/'manifest.json',manifest)


def test_portable_neighbor_cycle_closes_both_orders_at_all_cutoffs(package):
    result=module.replay(package)
    assert result['states']==14 and result['balanced_reactions']==14
    assert result['by_cutoff']['50']['deltaDeltaG_kJ_mol']==pytest.approx(-1.8616197516972124,abs=1e-8)
    assert not result['conformer_convergence_established']


@pytest.mark.parametrize('change',['scope','nan_G','duplicate','energy','cutoff','path_escape'])
def test_corrupt_or_misbound_neighbor_records_fail_even_if_rehashed(package,change):
    path=package/'manifest.json';manifest=json.loads(path.read_text())
    if change=='scope':manifest['electronic_coupling_computed']=True;write(path,manifest)
    elif change=='nan_G':manifest['states'][0]['checked_G_eV']=float('nan');write(path,manifest)
    elif change=='duplicate':manifest['states'][-1]['id']=manifest['states'][0]['id'];write(path,manifest)
    elif change=='energy':
        state=next(r for r in manifest['states'] if r['id']=='polymer_pair')
        rehash(package,state['files']['energy'],lambda d:d.update(e_smd_eV=d['e_smd_eV']+.1))
    elif change=='cutoff':
        rehash(package,manifest['neighbor_source'],lambda d:d['by_cutoff']['25'].update(deltaDeltaG_kJ_mol=float('nan')))
    else:
        manifest['states'][0]['files']['geometry']='../outside.xyz';write(path,manifest)
    with pytest.raises(ValueError):module.replay(package)
