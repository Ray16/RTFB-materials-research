from copy import deepcopy
import json
from pathlib import Path
import pytest
from redox.core.common import HARTREE_EV,EV_KJ
from redox.validation.component_sensitivity import reconstruct


def inputs():
    root=Path(__file__).resolve().parents[1]
    package=json.loads((root/'examples/exact_viologen_replay/data/manifest.json').read_text())
    rows={r['id']:dict(stored_Ha=-100.-i,installed_Ha=-100.-i,corrected_Ha=-100.-i,
                     charge=r['charge'],mult=r['mult']) for i,r in enumerate(package['states'])}
    return package,rows


def test_only_electronic_contribution_is_reported_with_correct_units_and_signs():
    package,rows=inputs();rows['graft_ox2']['corrected_Ha']+=1e-6
    result=reconstruct(package,rows);cycles={r['id']:r for r in result['cycles']}
    assert cycles['E1_V_vs_Fc']['electronic_contribution_to_shift']['corrected_minus_installed']==pytest.approx(HARTREE_EV*1e-6,abs=1e-12)
    assert cycles['graft_precursor_2_kJ_mol']['electronic_contribution_to_shift']['corrected_minus_installed']==pytest.approx(HARTREE_EV*EV_KJ*1e-6,abs=1e-10)
    assert result['total_free_energy_requalified'] is False and result['nuclear_terms_recomputed'] is False


@pytest.mark.parametrize('bad',['missing','wrong_spin','nonfinite','unbalanced_cycle'])
def test_incomplete_or_changed_physical_cycle_rejected(bad):
    package,rows=inputs()
    if bad=='missing':rows.pop('fc_neu')
    if bad=='wrong_spin':rows['graft_ox1']['mult']=1
    if bad=='nonfinite':rows['local']['installed_Ha']=float('nan')
    if bad=='unbalanced_cycle':package['reactions'][0]['coefficients']['graft_ox2']=-2
    with pytest.raises(ValueError):reconstruct(package,rows)
