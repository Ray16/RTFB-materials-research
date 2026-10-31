import importlib.util
from pathlib import Path

import pytest


@pytest.fixture
def closure():
    path=Path(__file__).resolve().parents[1]/'scripts/release/lock_xtb_runtime.py'
    spec=importlib.util.spec_from_file_location('tested_xtb_artifact_lock',path);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module.dependency_closure


def record(name,depends):
    return dict(name=name,version='1',build='build',url=f'https://conda.anaconda.org/conda-forge/linux-64/{name}.conda',
                md5='a'*32,sha256='b'*64,depends=depends)


def test_runtime_closure_includes_transitive_dependencies_and_stops_cycles(closure):
    records=dict(xtb=record('xtb',['fortran >=1','__glibc >=2.17']),fortran=record('fortran',['blas >=1']),blas=record('blas',['fortran >=1']),unrelated=record('unrelated',[]))
    assert [row['name'] for row in closure(records)]==['blas','fortran','xtb']


@pytest.mark.parametrize('change',['missing_dependency','private_URL','missing_SHA256'])
def test_incomplete_or_nonpublic_artifact_lock_is_refused(closure,change):
    records=dict(xtb=record('xtb',[]))
    if change=='missing_dependency':records['xtb']['depends']=['missing >=1']
    if change=='private_URL':records['xtb']['url']='https://user:secret@example.invalid/file.conda'
    if change=='missing_SHA256':records['xtb']['sha256']=''
    with pytest.raises(ValueError):closure(records)
