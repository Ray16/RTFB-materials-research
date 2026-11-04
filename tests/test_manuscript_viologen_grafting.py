import importlib.util
import json
from pathlib import Path
import shutil

import pytest

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('frozen_graft_checker',ROOT/'manuscript/check_viologen_grafting.py')
checker=importlib.util.module_from_spec(spec);spec.loader.exec_module(checker)


@pytest.fixture
def frozen(tmp_path):
    source=ROOT/'validation/viologen_grafting_2026-10-31'
    for a,b in [('reaction.json','viologen_grafting_result.json'),('values.tex','viologen_grafting_values.tex'),('QC_table.tex','viologen_grafting_QC_table.tex')]:
        shutil.copyfile(source/a,tmp_path/b)
    return tmp_path


def test_exact_graft_reconstructs_from_electronic_nuclear_and_spin_terms(frozen):
    result=checker.audit(frozen)
    assert result['deltaG_kJ_mol_by_cutoff']['50']==pytest.approx(-11.49332965)
    assert not result['material_admission']


@pytest.mark.parametrize('change',['energy','charge','base','promotion','missing','selection','macro','QC'])
def test_changed_chemistry_numerics_or_scope_are_rejected(frozen,change):
    p=frozen/'viologen_grafting_result.json';d=json.loads(p.read_text())
    if change=='energy':d['components']['product']['active_E_smd_eV']+=.01
    if change=='charge':d['components']['product']['charge']=1
    if change=='base':d['base_term_kJ_mol']=-70.
    if change=='promotion':d['material_admission']=True
    if change=='missing':d['precursor_conformers'].pop()
    if change=='selection':d['selected_precursor_conformer_index']=0
    if change=='QC':d['precursor_conformers'][0]['maximum_sampled_column_error']=.02
    if change=='macro':(frozen/'viologen_grafting_values.tex').write_text('changed')
    p.write_text(json.dumps(d))
    with pytest.raises(ValueError):checker.audit(frozen)
