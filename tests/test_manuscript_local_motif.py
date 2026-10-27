import copy
import importlib.util
import json
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('paper_claims', ROOT/'manuscript/check_claims.py')
claims = importlib.util.module_from_spec(spec)
spec.loader.exec_module(claims)
PANEL = json.loads((ROOT/'results/reviews/2026-10-26/local_motif_complete/summary.json').read_text())


def audit(panel):
    return claims.audit_local_motif(lambda name: panel)


def test_twelve_pair_energies_reconstruct():
    result = audit(PANEL)
    assert result['triplets'] == 12
    assert result['central_mECN_kJ_mol']['occ_2_3__pose0'] == pytest.approx(-24.76152457)
    assert not result['free_energy_computed']


@pytest.mark.parametrize('mutation', ['energy', 'duplicate', 'free_energy', 'basis'])
def test_rejects_changed_or_overstated_evidence(mutation):
    panel = copy.deepcopy(PANEL)
    if mutation == 'energy':
        panel['states'][0]['exact_state_energies_eV']['AB']['acetonitrile'] += 0.1
    elif mutation == 'duplicate':
        panel['states'][1] = copy.deepcopy(panel['states'][0])
    elif mutation == 'free_energy':
        panel['free_energies_computed'] = True
    else:
        panel['basis_superposition_correction_computed'] = True
    with pytest.raises(ValueError):
        audit(panel)
