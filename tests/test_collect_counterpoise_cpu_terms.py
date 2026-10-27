import importlib.util
import json
from pathlib import Path

import pytest


ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('cpu_cp_collect',ROOT/'scripts/validation/collect_counterpoise_cpu_terms.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)


def test_balanced_five_term_sum_and_scope_gates():
    core=module.api.load_module('cpu_cp_collect_test_protocol',ROOT/'scripts/validation/local_counterpoise.py')
    energy={'AB':-2.03,'A_own':-1.,'B_own':-1.,'A_ghost':-1.001,'B_ghost':-1.002}
    grid={'grid_coordinates_AU.npy':'a','grid_weights.npy':'b'}
    rows=[dict(role=role,E_Ha=value,settings=core.SETTINGS,common_grid_hashes=grid,stationarity_passed=True)
          for role,value in energy.items()]
    result=module.assemble(rows,core)
    assert result['counterpoise_gas_kJ_mol']==pytest.approx(-.027*core.HA_EV*core.EV_KJ)
    assert not result['solvent_counterpoise_computed'] and not result['thermal_free_energy_computed']
    with pytest.raises(ValueError,match='five unique'):
        module.assemble(rows[:-1]+[rows[0]],core)
    altered=[dict(row) for row in rows];altered[-1]['stationarity_passed']=False
    with pytest.raises(ValueError,match='stationarity'):
        module.assemble(altered,core)
    altered=[dict(row) for row in rows];altered[-1]['common_grid_hashes']={'grid_coordinates_AU.npy':'wrong','grid_weights.npy':'b'}
    with pytest.raises(ValueError,match='grid'):
        module.assemble(altered,core)


def test_five_saved_water_terms_reconstruct_reviewed_counterpoise():
    core=module.api.load_module('cpu_cp_collect_water_protocol',ROOT/'scripts/validation/local_counterpoise.py')
    base=ROOT/'calcs/validation/2026-10-26_counterpoise/CPU_fixed_auxiliary_factory'
    summary=json.loads((base/'summary.json').read_text());review=json.loads((ROOT/'results/reviews/2026-10-26/counterpoise_control/independent_review.json').read_text())
    rows=[dict(role=role,E_Ha=summary['terms'][role]['E_Ha'],settings=core.SETTINGS,
               common_grid_hashes=summary['common_grid_hashes'],stationarity_passed=True) for role in module.ROLES]
    result=module.assemble(rows,core)
    assert result['counterpoise_gas_kJ_mol']==pytest.approx(review['CPU']['counterpoise_gas_kJ_mol'],abs=1e-9)
