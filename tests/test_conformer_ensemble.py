import math

import numpy as np
import pytest

from redox.properties.conformer_ensemble import (
    EV_KJ_MOL, KB_EV_K, batch_change, cluster_minima, neighbor_cycle, partition, summarize,
)


def row(key, energy=-100, free=-99):
    return dict(key=key, E_eV=energy, G_eV=free)


def test_two_equally_populated_wells_have_RT_log_two_entropy():
    result = partition([-100000, -100000])
    assert result["G_eV"] == pytest.approx(-100000 - KB_EV_K * 298.15 * math.log(2))
    assert result["configurational_lowering_kJ_mol"] == pytest.approx(1.7182820758)
    assert result["conditional_populations"] == pytest.approx([0.5, 0.5])


def test_partition_is_shift_invariant_and_stable_for_large_energy_gaps():
    first = partition([0, 1, 100000])
    shifted = partition([-10000, -9999, 90000])
    assert shifted["G_eV"] - first["G_eV"] == pytest.approx(-10000)
    assert first["conditional_populations"] == shifted["conditional_populations"]
    assert sum(first["conditional_populations"]) == pytest.approx(1)


@pytest.mark.parametrize("bad", [None, True, float("nan"), float("inf")])
def test_uncomputed_or_nonfinite_free_energy_cannot_be_partitioned(bad):
    with pytest.raises(ValueError):
        partition([-100, bad])


def test_duplicate_search_hits_cannot_manufacture_configurational_entropy():
    records = [row("first"), row("repeat")]
    result = summarize(records, [["first", "repeat"]], ["first", "repeat"])
    assert result["G_eV"] == -99
    assert result["unique_minima"] == 1
    assert result["search_hits_are_weights"] is False


def test_missing_native_G_leaves_ensemble_blank_even_with_electronic_energies():
    result = summarize([row("a"), row("b", free=None)], [["a"], ["b"]], ["a", "b"])
    assert result["G_eV"] is None
    assert result["missing_keys"] == ["b"]


def test_groups_must_cover_each_expected_minimum_once():
    with pytest.raises(ValueError, match="partition"):
        summarize([row("a"), row("b")], [["a", "a"]], ["a", "b"])


def test_complete_link_clustering_prevents_transitive_chain_merges():
    records = [row("a"), row("b"), row("c")]
    distances = [[0, 0.08, 0.16], [0.08, 0, 0.08], [0.16, 0.08, 0]]
    assert cluster_minima(records, distances) == [["a", "b"], ["c"]]


def test_same_geometry_with_distinct_electronic_energy_is_not_merged():
    assert cluster_minima([row("a"), row("b", energy=-100 + 1 / EV_KJ_MOL)], np.zeros((2, 2))) == [["a"], ["b"]]


def test_batch_convergence_uses_added_wells_not_just_unchanged_minimum():
    records = [row("a"), row("b")]
    first = summarize(records[:1], [["a"]], ["a"])
    second = summarize(records, [["a"], ["b"]], ["a", "b"])
    result = batch_change(first, second, tolerance_kJ_mol=1)
    assert result["minimum_G_change_kJ_mol"] == 0
    assert result["ensemble_G_change_kJ_mol"] == pytest.approx(-1.7182820758)
    assert result["within_batch_tolerance"] is False
    assert result["global_ensemble_convergence_proven"] is False


def test_convergence_requires_strictly_nested_batches():
    first = summarize([row("a")], [["a"]], ["a"])
    with pytest.raises(ValueError, match="nested"):
        batch_change(first, first)


def test_neighbor_cycle_cancels_common_and_site_additive_terms():
    first = neighbor_cycle(-100, -110, -112, -121)
    shifted = neighbor_cycle(-10100, -10113, -10117, -10129)
    assert first == pytest.approx(EV_KJ_MOL)
    assert shifted == pytest.approx(first)
