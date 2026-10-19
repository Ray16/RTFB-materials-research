"""Spin-helper repairs must carry projected thermochemistry and roll back failed installs."""
import json
from pathlib import Path
import shutil

import pytest

from redox.qm import native_thermal_repair as repair
from redox.validation import spins


@pytest.fixture
def helper_repair(tmp_path,monkeypatch):
    helper,selected,review=[tmp_path/name for name in ('helper','selected','review')]
    for directory in (helper,selected,review):directory.mkdir()
    for directory,tag in ((helper,'old'),(selected,'new')):
        (directory/'opt.xyz').write_text(tag+' geometry')
        (directory/'result.json').write_text(json.dumps(dict(charge=0,mult=3)))
        for cache in ('sp','native_thermal','projected_thermal'):
            (directory/cache).mkdir();(directory/cache/f'{tag}.json').write_text(tag)
    (helper/'prior_failure.txt').write_text('preserve')
    before={str(p.relative_to(helper)):p.read_bytes() for p in helper.rglob('*') if p.is_file()}
    monkeypatch.setattr(repair,'select_policy',lambda directory:None)
    def qualified(directory,q,m):
        return all((Path(directory)/cache/'new.json').is_file()
                   for cache in ('sp','native_thermal','projected_thermal'))
    monkeypatch.setattr(repair,'load_record',qualified)
    monkeypatch.setattr(repair,'load_native_thermal',qualified)
    return helper,selected,review,before


def test_helper_install_preserves_complete_partition_and_old_evidence(helper_repair):
    helper,selected,review,before=helper_repair
    archive=spins._install_helper_minimum(helper,selected,review,0,3,spins.sha(helper/'opt.xyz'))
    assert {str(p.relative_to(archive)):p.read_bytes() for p in archive.rglob('*') if p.is_file()}==before
    assert (helper/'projected_thermal/new.json').read_text()=='new'
    assert (helper/'projected_thermal/old.json').read_text()=='old'
    assert (helper/'prior_failure.txt').read_text()=='preserve'
    assert (review/'install_producer_source.py').read_bytes()==repair.REPAIR_SOURCE_BYTES
    metadata=json.loads((helper/'result.json').read_text())
    assert metadata['minimum_repair_source']==str(archive/'opt.xyz')
    assert metadata['opt_out']==str(helper/'opt.xyz')


@pytest.mark.parametrize('failure',['missing_partition','stage_validation','post_install_validation'])
def test_failed_helper_install_leaves_original_bytes_intact(helper_repair,monkeypatch,failure):
    helper,selected,review,before=helper_repair
    if failure=='missing_partition':
        shutil.rmtree(selected/'projected_thermal')
    else:
        monkeypatch.setattr(repair,'load_native_thermal',lambda directory,q,m:
                            False if failure=='stage_validation' or Path(directory)==helper else True)
    with pytest.raises(RuntimeError):
        spins._install_helper_minimum(helper,selected,review,0,3,spins.sha(helper/'opt.xyz'))
    assert {str(p.relative_to(helper)):p.read_bytes() for p in helper.rglob('*') if p.is_file()}==before
    if failure=='post_install_validation':assert (review/'failed_install_preserved').is_dir()


def test_helper_source_change_during_refinement_blocks_installation(helper_repair):
    helper,selected,review,_=helper_repair
    source=spins.sha(helper/'opt.xyz')
    (helper/'opt.xyz').write_text('external replacement')
    current={str(p.relative_to(helper)):p.read_bytes() for p in helper.rglob('*') if p.is_file()}
    with pytest.raises(RuntimeError,match='source changed'):
        spins._install_helper_minimum(helper,selected,review,0,3,source)
    assert {str(p.relative_to(helper)):p.read_bytes() for p in helper.rglob('*') if p.is_file()}==current
