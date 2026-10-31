import importlib.util
import json
from pathlib import Path
import shutil

import pytest

ROOT=Path(__file__).resolve().parents[1]
definition=importlib.util.spec_from_file_location('checked_viologen_paper',ROOT/'manuscript/check_viologen_native.py')
checker=importlib.util.module_from_spec(definition);definition.loader.exec_module(checker)


@pytest.fixture
def frozen(tmp_path):
    source=ROOT/'validation/viologen_native_recovery_2026-10-30'
    for a,b in [('candidate.json','viologen_native_candidate.json'),('numeric_inputs.json','viologen_native_numeric_inputs.json'),
                ('native_QC.json','viologen_native_QC.json'),('figure.json','viologen_native_figure.json'),
                ('values.tex','viologen_native_values.tex'),('registry_identity_index.json','current_registry_identity_index.json'),
                ('figure.pdf','viologen_native_redox.pdf')]:shutil.copyfile(source/a,tmp_path/b)
    return tmp_path


def test_five_state_reconstruction_and_exact_registry_identity(frozen):
    result=checker.audit(frozen)
    assert result['potentials_V_vs_Fc']==pytest.approx([-.6077915087,-1.2080016581])
    assert result['registry_graphs']==904 and not result['material_admission']


@pytest.mark.parametrize('change',['energy','nonfinite','identity','basis','reference','QC','promotion','figure','macro','registry'])
def test_changed_components_scope_or_plot_cannot_enter_the_paper(frozen,change):
    p=frozen/'viologen_native_numeric_inputs.json';d=json.loads(p.read_text())
    if change=='energy':d['states'][0]['active_E_smd_eV']+=.1
    if change=='nonfinite':d['states'][0]['spin_partition_eV']=float('nan')
    if change=='identity':d['states'][0]['product_inchikey']='another'
    if change=='basis':d['states'][0]['active_protocol']['basis']='def2-tzvp'
    if change=='reference':d['states'][-1]['G_by_cutoff_eV']['50']+=.1
    if change=='QC':p=frozen/'viologen_native_QC.json';d=json.loads(p.read_text());d['states'][0]['final_minimum_frequency_cm']=-10.
    if change=='promotion':p=frozen/'viologen_native_candidate.json';d=json.loads(p.read_text());d['material_admission']=True
    if change=='registry':p=frozen/'current_registry_identity_index.json';d=json.loads(p.read_text());d['graphs']=d['graphs'][:-1]
    if change=='figure':(frozen/'viologen_native_redox.pdf').write_bytes(b'changed')
    if change=='macro':
        p2=frozen/'viologen_native_values.tex';p2.write_text(p2.read_text().replace('-0.608','-0.600'))
    p.write_text(json.dumps(d))
    with pytest.raises(ValueError):checker.audit(frozen)
