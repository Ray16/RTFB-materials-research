import importlib.util
import json
from pathlib import Path
import shutil

import pytest

ROOT=Path(__file__).resolve().parents[1]
definition=importlib.util.spec_from_file_location('frozen_tempo_solution',ROOT/'manuscript/check_tempo_solution.py')
checker=importlib.util.module_from_spec(definition);definition.loader.exec_module(checker)


@pytest.fixture
def frozen(tmp_path):
    source=ROOT/'validation/tempo_solution_pair_2026-10-30'
    for a,b in [('numeric_inputs.json','tempo_solution_numeric_inputs.json'),('independent_review.json','tempo_solution_review.json'),
                ('candidate.json','tempo_solution_candidate.json'),('figure.json','tempo_solution_figure.json'),
                ('figure.pdf','tempo_solution_relaxation.pdf'),
                ('values.tex','tempo_solution_values.tex')]:shutil.copyfile(source/a,tmp_path/b)
    return tmp_path


def test_six_cell_frozen_result_reconstructs_both_directions(frozen):
    result=checker.audit(frozen)
    assert result['reduction_eV']==pytest.approx(1.47712203)
    assert result['oxidation_eV']==pytest.approx(1.63864107)
    assert not result['symmetric_scalar_qualified']


@pytest.mark.parametrize('mutation',['energy','reference','role','spin','promotion','nonfinite','figure','macro'])
def test_changed_numeric_inputs_or_inference_cannot_pass(frozen,mutation):
    path=frozen/'tempo_solution_numeric_inputs.json';data=json.loads(path.read_text())
    if mutation=='energy':data['cells'][0]['e_eV']+=.1
    if mutation=='reference':
        row=next(r for r in data['cells'] if r['reference'] is not None);row['reference_density_sha256']='wrong'
    if mutation=='role':data['cells'][0]['at']='another geometry'
    if mutation=='spin':data['cells'][0]['mult']=3
    if mutation=='nonfinite':data['cells'][0]['e_Ha']=float('nan')
    if mutation=='promotion':
        path=frozen/'tempo_solution_candidate.json';data=json.loads(path.read_text());data['qualified_symmetric_single_body_eV']=1.558
    if mutation=='figure':(frozen/'tempo_solution_relaxation.pdf').write_bytes(b'changed')
    if mutation=='macro':
        macro=frozen/'tempo_solution_values.tex'
        macro.write_text(macro.read_text().replace('{1.477}','{1.478}'))
    path.write_text(json.dumps(data))
    with pytest.raises(ValueError):checker.audit(frozen)
