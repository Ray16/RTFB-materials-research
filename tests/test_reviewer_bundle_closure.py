import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess

import pytest

ROOT=Path(__file__).parents[1]
def load(name,path):
    spec=importlib.util.spec_from_file_location(name,ROOT/path);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module
prepare=load('bundle_v2_prepare','scripts/release/prepare_manuscript_bundle_v2.py')
validate=load('bundle_v2_validate','scripts/release/validate_manuscript_bundle_v2.py')


def make_project(root):
    for name in ['library/manifest.csv','src/redox/qm/xtb_runtime.py','scripts/pipeline/run_portable_thermal.py','results/control.json']:
        path=root/name;path.parent.mkdir(parents=True,exist_ok=True);path.write_text(name)
    control=root/'results/control.json';data={'evidence_sha256':{'results/control.json':hashlib.sha256(control.read_bytes()).hexdigest()}}
    data['certificate_content_sha256']=hashlib.sha256(json.dumps(data,sort_keys=True,separators=(',',':')).encode()).hexdigest()
    path=root/'config/projected_thermal_validation.json';path.parent.mkdir(parents=True);path.write_text(json.dumps(data))


def test_relocated_certificate_closure_refuses_missing_or_changed_dependency(tmp_path):
    make_project(tmp_path);result=validate.verify_closure(tmp_path)
    assert result['projected_certificate_dependencies']==1
    (tmp_path/'results/control.json').write_text('changed')
    with pytest.raises(ValueError,match='dependency differs'):validate.verify_closure(tmp_path)


def test_git_bundle_requires_all_certificate_sources_to_be_declared(tmp_path):
    make_project(tmp_path);subprocess.run(['git','init','-q',str(tmp_path)],check=True)
    subprocess.run(['git','add','.'],cwd=tmp_path,check=True)
    subprocess.run(['git','-c','user.name=Fixture','-c','user.email=fixture@example.invalid','commit','-qm','fixture'],cwd=tmp_path,check=True)
    head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=tmp_path,text=True).strip()
    config={'source_paths':['library','src','scripts','config','results']}
    assert prepare.check_git_closure(tmp_path,config,head)==1
    config['source_paths'].remove('results')
    with pytest.raises(ValueError,match='Missing projected certificate'):prepare.check_git_closure(tmp_path,config,head)
    assert not prepare.included('library_extra/manifest.csv',['library'])
    assert not prepare.included('../library/manifest.csv',['library'])
