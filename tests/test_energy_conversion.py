from decimal import Decimal
import json
import math
from pathlib import Path

import pytest

from redox.core.common import EV_KJ, LEGACY_EV_KJ
from redox.properties.conformer_ensemble import EV_KJ_MOL
from redox.workflows.polymer_campaign import conditional_grafting


def test_one_ev_per_molecule_equals_the_exact_SI_molar_conversion():
    # Independent decimal arithmetic from the two exact SI defining values.
    expected = Decimal('1.602176634e-19')*Decimal('6.02214076e23')/Decimal('1000')
    assert EV_KJ == float(expected)
    assert EV_KJ == EV_KJ_MOL
    assert LEGACY_EV_KJ == 96.485 and EV_KJ != LEGACY_EV_KJ


def test_cached_reaction_sums_use_exact_units_without_rewriting_the_source():
    path = Path(__file__).resolve().parents[1]/'examples/cached_thermochemistry/data/manifest.json'
    source_bytes = path.read_bytes()
    data = json.loads(source_bytes)
    states = {r['id']: r['checked_G_eV'] for r in data['states']}
    expected = {'HCl_reference': -0.2915709641, 'DBU_separated_ions': -103.0688261726}
    for reaction in data['reactions']:
        delta = math.fsum(states[key]*coefficient for key, coefficient in reaction['coefficients'].items())
        assert delta*EV_KJ == pytest.approx(expected[reaction['name']], abs=1e-7)
        assert delta*LEGACY_EV_KJ == pytest.approx(reaction['source_kJ_mol'], abs=1e-7)
        assert abs(delta*(EV_KJ-LEGACY_EV_KJ)) < .0004
    assert path.read_bytes() == source_bytes


def test_conditional_cycle_reports_its_conversion_and_cancels_shared_terms():
    result = conditional_grafting(-100., -110., -112., -121., 5.)
    shifted = conditional_grafting(-100., -110., -112., -121., 8.)
    assert result['conversion_kJ_mol_per_eV'] == EV_KJ
    assert result['neighbor_penalty_deltaDeltaG_kJ_mol'] == pytest.approx(EV_KJ)
    assert shifted['neighbor_penalty_deltaDeltaG_kJ_mol'] == result['neighbor_penalty_deltaDeltaG_kJ_mol']
    assert (result['graft_j_given_i_kJ_mol']-result['graft_j_on_empty_kJ_mol']) == pytest.approx(EV_KJ)
