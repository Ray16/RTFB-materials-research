import importlib.util
from pathlib import Path
from types import SimpleNamespace

import pytest

ROOT=Path(__file__).parents[1]
spec=importlib.util.spec_from_file_location('CPU_dependency_test',ROOT/'scripts/fleet/preflight_cpu_dependencies.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)


@pytest.mark.parametrize('name,distribution',[('gpu4pyscf','gpu4pyscf-cuda12x'),('cupy','cupy-cuda12x'),
    ('torch','torch'),('fairchem.core','fairchem-core')])
def test_GPU_dependencies_are_never_imported_in_CPU_preflight(name,distribution):
    calls=[]
    def forbidden(name):raise AssertionError('A CPU probe attempted a GPU import')
    def version(name):calls.append(name);return 'test-version'
    assert module.module_version(name,forbidden,version)==('test-version',True)
    assert calls==[distribution]


def test_CPU_dependency_is_actually_imported():
    calls=[]
    def importer(name):calls.append(name);return SimpleNamespace(__version__='checked')
    assert module.module_version('numpy',importer,lambda name:None)==('checked',False)
    assert calls==['numpy']


def test_missing_GPU_distribution_is_not_silently_accepted():
    def missing(name):raise module.metadata.PackageNotFoundError(name)
    with pytest.raises(module.metadata.PackageNotFoundError):module.module_version('gpu4pyscf',version=missing)
