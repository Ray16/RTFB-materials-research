import importlib.util
import inspect
from pathlib import Path

import numpy as np
import pytest

ROOT=Path(__file__).parents[1]
spec=importlib.util.spec_from_file_location('smooth_geometry_schema_test',ROOT/'scripts/pipeline/run_smooth_geometry.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)


def test_geometry_adapter_cannot_emit_primary_schema(tmp_path):
    body="""def run(sd):
    forcecheck=forces.final_energy_check(sd/'opt.xyz',sd/'force')
    settings=dict(convergence_gmax=1e-4,convergence_grms=6e-5,convergence_energy=1e-7,
                        convergence_dmax=.002,convergence_drms=.001)
    atomic_write_json(sd/'result.json',settings)
    return sd/'result.json',forcecheck
"""
    writes=[];local={};scope=dict(atomic_write_json=lambda path,row:writes.append((path,row)),
        forces=None,_expanded_final_check=lambda *args:'independent force check')
    exec(compile(module.alternate_run_source(body),'<schema control>','exec'),scope,local)
    path,check=local['run'](tmp_path)
    assert path.name=='experimental_geometry_result.json'
    assert writes[0][0].name=='experimental_geometry_result.json'
    assert writes[0][1]['convergence_gmax']==1e-5
    assert check=='independent force check'


@pytest.mark.parametrize('axes',[[],[(0,1),(1,2)]])
def test_critical_force_checks_preserve_largest_components_and_integer_indices(axes):
    text='''def selectors(force):
    indices = np.argsort(abs(force).ravel())[-2:]
    return indices
'''
    namespace=dict(np=np);exec(module.expanded_force_source(text,axes),namespace)
    indices=namespace['selectors'](np.arange(6.).reshape(2,3))
    assert indices.dtype.kind in 'iu'
    assert {4,5}.issubset(indices)
    assert {3*a+b for a,b in axes}.issubset(indices)


def test_unknown_geometry_writer_template_requires_review():
    with pytest.raises(ValueError,match='primary-schema'):
        module.alternate_run_source('def run(): pass')
