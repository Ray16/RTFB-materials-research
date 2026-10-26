import copy
import importlib.util
import json
from pathlib import Path

import pytest


ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('manuscript_tempo_audit',ROOT/'manuscript/check_claims.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)


def evidence():
    root=ROOT/'manuscript/evidence/2026-10-26_tempo_exact_support'
    names=['tempo_route.json','tempo_route_screen.json','tempo_native_potential.json','tempo_decision.json','tempo_summary.json']
    return {name:json.loads((root/name).read_text()) for name in names}


def test_tempo_claims_reconstruct_and_remain_pending():
    data=evidence();r=module.audit_tempo_exact_support(lambda name:copy.deepcopy(data[name]))
    assert r['route']=='proceed' and r['material_decision']=='pending'
    assert r['potential_V_vs_Fc_by_cutoff']['50']==pytest.approx(.5787750407616841)


@pytest.mark.parametrize('mutation',['potential','route_energy','eligible'])
def test_tempo_claim_audit_rejects_misstatements(mutation):
    data=evidence()
    if mutation=='potential':data['tempo_native_potential.json']['potential_V_vs_Fc_by_rotor_cutoff']['50']+=.1
    if mutation=='route_energy':data['tempo_summary.json']['actual_route_deltaG_computed']=True
    if mutation=='eligible':data['tempo_decision.json']['decision']['all_required_axes_pass']=True
    with pytest.raises(ValueError):module.audit_tempo_exact_support(lambda name:copy.deepcopy(data[name]))
