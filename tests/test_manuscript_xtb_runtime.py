import copy
import importlib.util
import json
from pathlib import Path

import pytest


@pytest.fixture
def evidence(tmp_path):
    root=Path(__file__).resolve().parents[1]
    spec=importlib.util.spec_from_file_location('checked_manuscript_xtb',root/'manuscript/check_xtb_runtime.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    source=root/'validation/xtb_runtime_relocation_2026-10-31'
    rows={target:json.loads((source/original).read_text()) for original,target in [('summary.json','xtb_runtime_installation.json'),('packages.json','xtb_runtime_packages.json'),('installed_packages.json','xtb_installed_packages.json')]}
    return module,rows,tmp_path


def write(rows,directory):
    for name,data in rows.items():(directory/name).write_text(json.dumps(data))


def test_frozen_component_check_reconstructs_the_methods_scope(evidence):
    module,rows,directory=evidence;write(rows,directory)
    assert module.audit(directory)['packages']==18


@pytest.mark.parametrize('change',['kernel','soft_mode','mask','build','scope','CUDA','artifact'])
def test_changed_scientific_or_installation_evidence_is_rejected(evidence,change):
    module,rows,directory=evidence;rows=copy.deepcopy(rows);run=rows['xtb_runtime_installation.json']
    if change=='kernel':run['numerical_checks']['upstream_kernel']['G_Ha']+=1e-5
    if change=='soft_mode':run['numerical_checks']['soft_mode_retained']=False
    if change=='mask':run['original_xtb_files_masked']=False
    if change=='build':next(r for r in rows['xtb_runtime_packages.json']['packages'] if r['name']=='xtb')['build']='other'
    if change=='scope':run['full_quantum_environment_qualified']=True
    if change=='CUDA':run['numerical_checks']['CUDA_imported']=True
    if change=='artifact':rows['xtb_installed_packages.json']['xtb']['sha256']='changed'
    write(rows,directory)
    with pytest.raises(ValueError):module.audit(directory)
