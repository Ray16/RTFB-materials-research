"""Frozen publication values must preserve state, translation and duplicate scope."""
import importlib.util
import json
from pathlib import Path
import shutil

import pytest

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('counterion_publication_audit',ROOT/'manuscript/check_counterion_association.py')
checker=importlib.util.module_from_spec(spec);spec.loader.exec_module(checker)


@pytest.fixture(params=['legacy','expanded'])
def frozen(tmp_path,request):
    source=ROOT/('validation/counterion_sensitivity_2026-11-02/publication' if request.param=='legacy'
                else 'validation/counterion_completion_2026-11-03/publication')
    target=tmp_path/'frozen';shutil.copytree(source,target)
    return target


def edit(folder,name,change):
    path=folder/name;data=json.loads(path.read_text());change(data);path.write_text(json.dumps(data))


def test_actual_frozen_outcomes_replay_with_one_statistical_group_per_near_duplicate(frozen):
    result=checker.audit(frozen)
    count=json.loads((frozen/'counterion_inputs.json').read_text())['checked_outcomes']
    assert result['checked_outcomes']==count and result['duplicate_groups']==count-1
    low,high=(15.41,15.43) if count==4 else (5.32,5.34)
    assert low<result['minimum_primary_association_1M_kJ_mol']<high
    assert 22.35<result['maximum_primary_association_1M_kJ_mol']<22.37
    assert result['same_group_N0_heavy_RMSD_A']<.004
    assert not result['conformer_ensemble_qualified'] and not result['experimental_accuracy_validated']


@pytest.mark.parametrize('change',['redox_spin','energy','thermal_definition','model','lost_chloride',
                                  'double_conversion','duplicate_weights','population_claim'])
def test_changed_chemistry_or_unsupported_interpretation_is_rejected(frozen,change):
    report='counterion_review_N0_minus.json'
    if change=='redox_spin':edit(frozen,report,lambda d:d.update(mult=2))
    if change=='energy':edit(frozen,report,lambda d:d['G_by_cutoff_eV'].update({'50':d['G_by_cutoff_eV']['50']+.01}))
    if change=='thermal_definition':edit(frozen,report,lambda d:d['thermal_protocol'].update(temperature_K=300))
    if change=='model':edit(frozen,report,lambda d:d['electronic_protocol'].update(solvent_model='other'))
    if change=='lost_chloride':edit(frozen,report,lambda d:d.update(smiles=d['smiles'].replace('.[Cl-]','')))
    if change=='double_conversion':
        def twice(d):
            row=next(r for r in d['states'] if r['key']=='N0_above')
            row['one_M_deltaG_kJ_mol']['50']+=row['one_atm_to_one_M_correction_kJ_mol']
        edit(frozen,'counterion_ledger_N0_minus.json',twice)
    if change=='duplicate_weights':
        edit(frozen,'counterion_duplicate_review.json',lambda d:d['groups_by_RMSD_cutoff_A'].update({'0.1':[['N0_minus'],['N0_plus'],['N1_above'],['N1_below']]}))
    if change=='population_claim':edit(frozen,'counterion_inputs.json',lambda d:d.update(conformer_ensemble_qualified=True))
    with pytest.raises(ValueError):checker.audit(frozen)
