from copy import deepcopy
import json
from pathlib import Path

import pytest

from redox.core.common import EV_KJ
from redox.validation.native_component_sensitivity import reconstruct


@pytest.fixture
def cohort():
    package_path=Path(__file__).resolve().parents[1]/'examples/exact_viologen_replay/data'
    manifest=json.loads((package_path/'manifest.json').read_text());rows={}
    context=dict(quantity='composite_solution_G',unit='eV',electronic_protocol=manifest['active_protocol'],
        thermal_definition=dict(native_protocol=manifest['native_protocol']),solvent='acetonitrile',solvent_model='SMD',
        temperature_K=298.15,standard_state='same source translation convention',reference='same declared primary model')
    for state in manifest['states']:
        lines=(package_path/state['files']['geometry']).read_text().splitlines();symbols=[line.split()[0] for line in lines[2:2+int(lines[0])]]
        rows[state['id']]=dict(charge=state['charge'],mult=state['mult'],qc='independently_checked',
            sources={'fixture':'arithmetic test only'},native_component_requalified=True,symbols=symbols,
            context=deepcopy(context),G_by_cutoff_eV=deepcopy(state['checked_G_by_cutoff_eV']))
    return manifest,rows


def test_complete_cohort_reconstructs_original_cycles_and_retains_scope(cohort):
    result=reconstruct(*cohort)
    assert result['states']==10 and len(result['cycles'])==6
    assert result['total_free_energy_cycles_numerically_requalified']
    assert not result['experimental_accuracy_validated'] and not result['material_admission']
    assert max(abs(v) for r in result['cycles'] for v in r['corrected_minus_source_by_cutoff'].values())<1e-9


def test_redox_sign_units_and_grafting_shift_follow_declared_stoichiometry(cohort):
    manifest,rows=cohort
    for c in ('25','50','100'):rows['graft_ox2']['G_by_cutoff_eV'][c]+=.001
    cycles={r['id']:r for r in reconstruct(manifest,rows)['cycles']}
    assert cycles['E1_V_vs_Fc']['corrected_minus_source_by_cutoff']['50']==pytest.approx(.001,abs=1e-10)
    assert cycles['graft_precursor_0_kJ_mol']['corrected_minus_source_by_cutoff']['50']==pytest.approx(.001*EV_KJ,abs=1e-8)


@pytest.mark.parametrize('bad',['missing','state','atoms','nonfinite','pending','unit','quantity','temperature','solvent','thermal','reference','coefficients','empty'])
def test_incomplete_mixed_or_mislabeled_G_cohort_is_rejected(cohort,bad):
    manifest,rows=cohort;r=rows['graft_ox2']
    if bad=='missing':rows.pop('fc_neu')
    if bad=='state':r['mult']=3
    if bad=='atoms':r['symbols'].append('H')
    if bad=='nonfinite':r['G_by_cutoff_eV']['50']=float('nan')
    if bad=='pending':r['native_component_requalified']=False
    if bad=='unit':r['context']['unit']='hartree'
    if bad=='quantity':r['context']['quantity']='electronic_energy'
    if bad=='temperature':r['context']['temperature_K']=310.
    if bad=='solvent':r['context']['solvent']='water'
    if bad=='thermal':r['context']['thermal_definition']['native_protocol']['basis']='different'
    if bad=='reference':r['context']['reference']='other'
    if bad=='coefficients':manifest['reactions'][0]['coefficients']['graft_ox2']=-.5
    if bad=='empty':manifest['reactions']=[]
    with pytest.raises(ValueError):reconstruct(manifest,rows)
