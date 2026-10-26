import importlib.util
import math
from pathlib import Path

import pytest


script = Path(__file__).resolve().parents[1]/'scripts/validation/empty_conformer_subset_sensitivity.py'
spec = importlib.util.spec_from_file_location('empty_subsets', script)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def test_two_equal_minima_have_rt_log_two_subset_effect():
    rows = [dict(key=k, G_by_cutoff_eV={'50': -10.0}) for k in ('a', 'b')]
    result = module.subset_diagnostics(rows, '50')
    group = result['subsets_by_size']['1']
    expected = 8.31446261815324 * 298.15 * math.log(2) / 1000
    assert group['combinations'] == 2
    assert group['max_partition_difference_kJ_mol'] == pytest.approx(expected, abs=1e-9)
    assert group['max_minimum_difference_kJ_mol'] == 0


def test_exhaustive_subsets_identify_missed_lowest_minimum():
    rows = [dict(key=k, G_by_cutoff_eV={'50': g}) for k, g in [('a', 0.0), ('b', 0.1), ('c', 0.2)]]
    result = module.subset_diagnostics(rows, '50')
    group = result['subsets_by_size']['2']
    assert group['combinations'] == 3
    worst = max(group['subsets'], key=lambda r: r['subset_minus_full_minimum_G_kJ_mol'])
    assert worst['keys'] == ['b', 'c']
    assert worst['subset_minus_full_minimum_G_kJ_mol'] == pytest.approx(0.1*module.EV_KJ_MOL)
    assert all(row['subset_minus_full_G_kJ_mol'] >= 0 for row in group['subsets'])


def test_rejects_duplicate_identity_and_nonfinite_free_energy():
    with pytest.raises(ValueError):
        module.subset_diagnostics([{'key':'x'}, {'key':'x'}], '50')
    with pytest.raises(ValueError):
        module.subset_diagnostics([dict(key='x', G_by_cutoff_eV={'50':0}),
                                   dict(key='y', G_by_cutoff_eV={'50':float('nan')})], '50')
