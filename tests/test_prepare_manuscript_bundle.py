import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess

import pytest

ROOT=Path(__file__).parents[1]
spec=importlib.util.spec_from_file_location('prepare_manuscript_bundle',ROOT/'scripts/release/prepare_manuscript_bundle.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)


@pytest.fixture
def repository(tmp_path):
    root=tmp_path/'repo';root.mkdir()
    subprocess.run(['git','init','-q',str(root)],check=True)
    (root/'source.py').write_text('print("checked")\n')
    archive=root/'a.tar.zst';archive.write_bytes(b'previously checked compressed archive')
    record=root/'record.json';record.write_text(json.dumps({'compressed_archive':'a.tar.zst','archive_sha256':hashlib.sha256(archive.read_bytes()).hexdigest(),'scope':'fixture'}))
    cfg=root/'config.json';cfg.write_text(json.dumps({'schema':1,'source_paths':['source.py'],
        'packages':[{'id':'one','record':'record.json'}], 'public_deposition_authorized_by_manifest':False,
        'data_license':None,'author_metadata':None,'third_party_source_policy':'separate'}))
    subprocess.run(['git','add','source.py'],cwd=root,check=True)
    subprocess.run(['git','-c','user.name=Fixture','-c','user.email=fixture@example.invalid','commit','-qm','fixture'],cwd=root,check=True)
    return root,cfg,archive


def test_bundle_binds_committed_source_and_checked_archive_without_claiming_publication(repository,tmp_path):
    root,cfg,archive=repository;out=tmp_path/'bundle';result=module.prepare(root,cfg,out)
    assert result['status']=='local_manuscript_reviewer_bundle_prepared_not_public'
    assert result['files_sha256']['archives/a.tar.zst']==hashlib.sha256(archive.read_bytes()).hexdigest()
    assert not result['public_archive_deposited'] and not result['experimental_accuracy_inferred']
    assert (out/'project-source.tar').is_file()
    with pytest.raises(FileExistsError):module.prepare(root,cfg,out)


def test_changed_archive_and_escaping_source_are_refused_before_output(repository,tmp_path):
    root,cfg,archive=repository;archive.write_bytes(b'changed');out=tmp_path/'bundle'
    with pytest.raises(ValueError,match='archive hash'):module.prepare(root,cfg,out)
    assert not out.exists()
    d=json.loads(cfg.read_text());d['source_paths']=['../outside'];cfg.write_text(json.dumps(d))
    with pytest.raises(ValueError,match='Unsafe'):module.prepare(root,cfg,out)
    assert not out.exists()
