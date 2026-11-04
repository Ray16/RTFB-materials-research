"""Equivalent source declarations are safe; held routes and conflicting states are not."""
from copy import deepcopy
import json

import pytest

from redox.core.common import ROOT
from redox.build.contracts import state_declaration
from redox.workflows.route_molecular_panel import route_candidate

POOL=ROOT/'results/candidates/pool/2026-11-02_ptz_precursor_access'


def baseline():
    data=json.loads((POOL/'registry.json').read_text())
    return next(r for r in data['candidates'] if 'phenothiazine' in r['aliases'])


def test_repeated_exact_baseline_declarations_retain_one_unambiguous_state_contract():
    candidate=baseline();rows=[r for r in candidate['declarations'] if r['id']=='phenothiazine']
    assert len(rows)==2 and rows[0]['fragment']!=rows[1]['fragment']
    declaration=state_declaration(candidate,'phenothiazine')
    assert declaration['states']==[['neu',0,1,0],['ox',1,2,1]]


@pytest.mark.parametrize('field,value',[('states',[['neu',0,1,0],['ox',1,4,1]]),
    ('precursor_smiles','CCO'),('fragment','[*:1]OC')])
def test_conflicting_charge_spin_or_chemistry_cannot_be_hidden_by_a_matching_first_source(field,value):
    candidate=baseline();candidate['declarations'][-1][field]=value
    with pytest.raises(ValueError):state_declaration(candidate,'phenothiazine')


def test_route_passes_only_for_the_exact_reviewed_baseline():
    candidate,declaration,gate=route_candidate(POOL/'registry.json',POOL/'inputs/route_evidence.json',
        POOL/'inputs/assessment.json','phenothiazine','phenothiazine')
    assert gate['proceed_to_full_redox'] and gate['exact_route_binding_verified'] and not candidate['evaluations']
    with pytest.raises(ValueError,match='route must pass'):
        route_candidate(POOL/'registry.json',POOL/'inputs/route_evidence.json',POOL/'inputs/assessment.json',
            'proposal_ptz_176c4d2f397734ea','proposal_ptz_176c4d2f397734ea')
