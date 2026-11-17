import importlib.util
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
specification = importlib.util.spec_from_file_location('native_batch_partition_policy',
    ROOT/'scripts/validation/analyze_native_conformer_batches.py')
analysis = importlib.util.module_from_spec(specification)
specification.loader.exec_module(analysis)


def row(index, free=-99., energy=-100.):
    return dict(key=f'all_same__occ_none__conf{index:03d}', E_eV=energy,
                G_by_cutoff_eV={str(rotor): free for rotor in (25, 50, 100)})


def test_added_duplicate_search_hit_cannot_lower_ensemble_free_energy():
    rows = [row(1), row(2, free=-98.9, energy=-99.9), row(3)]
    matrix = np.array([[0, 2, .04], [2, 0, 2], [.04, 2, 0]])
    result = analysis.finite_batches(rows, matrix, first_count=2)
    for sensitivity in result.values():
        for rotor in sensitivity.values():
            assert rotor['first']['unique_minima'] == rotor['expanded']['unique_minima'] == 2
            assert rotor['comparison']['ensemble_G_change_kJ_mol'] == pytest.approx(0)
            assert rotor['expanded']['search_hits_are_weights'] is False


def test_added_distinct_well_changes_partition_even_when_minimum_is_unchanged():
    rows = [row(1), row(2), row(3)]
    matrix = np.ones((3, 3))-np.eye(3)
    result = analysis.finite_batches(rows, matrix, first_count=2)['0.1']['50']
    assert result['comparison']['minimum_G_change_kJ_mol'] == 0
    assert result['comparison']['ensemble_G_change_kJ_mol'] < 0
    assert result['comparison']['global_ensemble_convergence_proven'] is False


def test_rmsd_sensitivity_is_an_explicit_scenario_not_an_error_bar():
    rows = [row(1), row(2), row(3)]
    matrix = np.array([[0, 2, .15], [2, 0, 2], [.15, 2, 0]])
    result = analysis.finite_batches(rows, matrix, first_count=2)
    assert result['0.1']['50']['expanded']['unique_minima'] == 3
    assert result['0.2']['50']['expanded']['unique_minima'] == 2
