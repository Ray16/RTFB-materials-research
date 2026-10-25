from copy import deepcopy
import importlib.util
import json
from pathlib import Path

import pytest

ROOT=Path(__file__).parents[1]
spec=importlib.util.spec_from_file_location('surface_claims_test',ROOT/'manuscript/check_claims.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)


@pytest.fixture
def data():
    folder=ROOT/'manuscript/evidence/2026-10-25_surface_controls'
    names=['CDS_contact.json','CDS_smooth_component.json','CDS_smooth_qualification.json',
        'CDS_smooth_full_gradients.json','CDS_smooth_quantum_CHO.json','CDS_smooth_quantum_NCl.json']
    return {name:json.loads((folder/name).read_text()) for name in names}


def test_completed_surface_control_scope_reconstructs(data):
    value=module.audit_surface_controls(data.__getitem__)
    assert value['neutral_small_quantum_controls']==4
    assert value['all_axis_component_states']==4
    assert value['polymer_G_qualification_inferred'] is False


@pytest.mark.parametrize('error',['energy','polymer_G','gradient_coverage','Hessian_coverage','model','nan'])
def test_changed_surface_evidence_cannot_support_paper_claims(data,error):
    data=deepcopy(data)
    if error=='energy':next(r for r in data['CDS_contact.json']['points'] if r['offset_A']==.001)['E_CDS_Ha']+=.01
    elif error=='polymer_G':data['CDS_smooth_qualification.json']['full_molecular_polymer_G_qualified']=True
    elif error=='gradient_coverage':data['CDS_smooth_full_gradients.json']['states']['base']['checks'].pop()
    elif error=='Hessian_coverage':data['CDS_smooth_quantum_CHO.json']['controls'][0]['arms']['smooth']['tested_Hessian_columns']=12
    elif error=='model':data['CDS_smooth_quantum_NCl.json']['portable_method_id']='another model'
    else:data['CDS_smooth_quantum_CHO.json']['controls'][0]['arms']['smooth']['maximum_force_error_Ha_Bohr']=float('nan')
    with pytest.raises(ValueError):module.audit_surface_controls(data.__getitem__)
