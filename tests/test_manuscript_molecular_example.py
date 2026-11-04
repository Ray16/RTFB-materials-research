import importlib.util
import json
from pathlib import Path
import shutil

import pytest

ROOT=Path(__file__).parents[1]
spec=importlib.util.spec_from_file_location('paper_molecular_example',ROOT/'manuscript/check_molecular_example.py')
checker=importlib.util.module_from_spec(spec);spec.loader.exec_module(checker)


@pytest.fixture
def frozen(tmp_path):
    folder=ROOT/'validation/exact_molecular_CPU_2026-10-31'
    for source,name in [('validation.json','exact_molecular_CPU_validation.json'),('replay.json','exact_molecular_CPU_replay.json'),('package_manifest.json','exact_molecular_CPU_manifest.json')]:shutil.copyfile(folder/source,tmp_path/name)
    return tmp_path


def test_versioned_example_agrees_with_extracted_fresh_environment(frozen):
    assert checker.audit(frozen,ROOT)['states']==10


@pytest.mark.parametrize('change',['scope','count','dependency','value','manifest'])
def test_changed_proof_scope_or_reported_values_cannot_enter_paper(frozen,change):
    p=frozen/'exact_molecular_CPU_validation.json';data=json.loads(p.read_text())
    if change=='scope':data['full_historical_source_closure_revalidated']=True
    if change=='count':data['states']=9
    if change=='dependency':data['scientific_dependency']='unlocked'
    if change=='value':
        p=frozen/'exact_molecular_CPU_replay.json';data=json.loads(p.read_text());data['by_cutoff']['50']['E1_V_vs_Fc']+=.01
    if change=='manifest':
        p=frozen/'exact_molecular_CPU_manifest.json';data=json.loads(p.read_text());data['material_admission']=True
    p.write_text(json.dumps(data))
    with pytest.raises(ValueError):checker.audit(frozen,ROOT)
