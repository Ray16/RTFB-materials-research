"""A publication benchmark cannot hide missing endpoints or create independence."""

import copy
import math

import pytest

from redox.core.common import load_config
from redox.validation.benchmark_panel import (
    cluster_intervals,
    freeze,
    key,
    metrics,
    number,
    summarize,
    validate_panel,
)


@pytest.fixture
def panel(tmp_path):
    source = tmp_path / "source.csv"
    source.write_text("frozen source")
    return freeze(load_config("benchmark_mecn").BENCHMARK, [source])


def comparison(target, residual=0.1):
    value = target["exp_V"] + residual
    published = round(value, 3)
    return dict(
        target,
        kind="E",
        calc_V=published,
        calc_unrounded_V=value,
        error_V=published - target["exp_V"],
        status="checked_native_parent",
    )


def test_entire_strict_catalog_is_frozen_without_outcome_selection(panel):
    assert len(panel["targets"]) == 39
    assert len({r["id"] for r in panel["targets"]}) == 30
    assert not panel["held_out_from_model_development"]
    assert all(
        r["tier"] == "A" and not r["excluded_molecule"] for r in panel["targets"]
    )
    assert any(r["tier"] != "A" for r in panel["source_catalog"])
    validate_panel(panel)


def test_changed_source_or_target_invalidates_panel(panel):
    changed = copy.deepcopy(panel)
    changed["targets"][0]["exp_V"] += 0.01
    with pytest.raises(ValueError, match="panel content changed"):
        validate_panel(changed)
    from pathlib import Path

    Path(next(iter(panel["source_sha256"]))).write_text("new source")
    with pytest.raises(ValueError, match="Changed benchmark evidence"):
        validate_panel(panel)


def test_missing_values_are_retained_in_coverage_and_never_zero_filled(panel):
    target = panel["targets"][0]
    report = summarize(panel, [comparison(target)])
    assert report["status"] == "incomplete_benchmark"
    assert report["admitted_events"] == 1 and report["pending_events"] == 38
    assert report["metrics"]["MAE_V"] == pytest.approx(0.1)
    assert sum(r["target_events"] for r in report["coverage_by_family"].values()) == 39
    assert all(u["intervals"] is None for u in report["uncertainty"])


def test_rounding_is_separate_from_primary_unrounded_residual(panel):
    row = comparison(panel["targets"][0], 0.1234567)
    report = summarize(panel, [row])
    assert report["metrics"]["MAE_V"] == pytest.approx(0.1234567)
    assert report["outliers"][0]["serialized_calc_V"] == row["calc_V"]


@pytest.mark.parametrize(
    "field,value",
    [
        ("tier", "B"),
        ("exp_V", 0.98765),
        ("error_V", 100),
        ("calc_unrounded_V", None),
        ("calc_unrounded_V", float("nan")),
        ("status", "optimizer exited zero"),
    ],
)
def test_changed_or_unqualified_endpoint_fails_closed(panel, field, value):
    row = comparison(panel["targets"][0])
    row[field] = value
    with pytest.raises(ValueError):
        summarize(panel, [row])


def test_duplicate_endpoint_cannot_inflate_sample_size(panel):
    row = comparison(panel["targets"][0])
    with pytest.raises(ValueError, match="Repeated comparison"):
        summarize(panel, [row, row])


def test_second_wave_is_not_renumbered_when_first_is_missing(panel):
    second = next(r for r in panel["targets"] if r["wave"] == 2)
    report = summarize(panel, [comparison(second)])
    assert set(report["by_wave"]) == {"2"}
    assert report["signed_wave_spacings"] == []


def test_signed_oxidation_spacing_is_reference_free_and_unrounded(panel):
    pair = [
        r
        for r in panel["targets"]
        if r["id"] == "bm_n_ethyl_3_7_dimethoxyphenothiazine"
    ]
    rows = [comparison(r, 0.1234567 if r["wave"] == 1 else 0.2345678) for r in pair]
    report = summarize(panel, rows)
    spacing = report["signed_wave_spacings"][0]
    assert spacing["exp_V"] == pytest.approx(-0.599)
    assert spacing["error_V"] == pytest.approx(-0.1111111)
    assert report["admitted_molecules"] == report["original_studies"] == 1
    assert report["molecule_macro_MAE_V"] == pytest.approx((0.1234567 + 0.2345678) / 2)


def test_two_waves_and_repeated_sources_are_not_independent_molecules(panel):
    pair = [
        r
        for r in panel["targets"]
        if r["id"] == "bm_n_ethyl_3_7_dimethoxyphenothiazine"
    ]
    other = next(r for r in panel["targets"] if r["id"] != pair[0]["id"])
    report = summarize(
        panel, [*(comparison(r, 0.1) for r in pair), comparison(other, 0.4)]
    )
    assert report["metrics"]["MAE_V"] == pytest.approx(0.2)
    assert report["molecule_macro_MAE_V"] == pytest.approx(0.25)


def test_cluster_draw_retains_all_members_and_is_deterministic():
    rows = [dict(id=str(i), error_V=i / 10) for i in range(5)]
    report = cluster_intervals(rows, "id", replicates=100)
    assert report == cluster_intervals(rows, "id", replicates=100)
    assert report["independent_cluster_count"] == 5
    assert report["intervals"]["MAE_V"][0] < 0.2 < report["intervals"]["MAE_V"][1]
    same = cluster_intervals(
        [dict(id=str(i), error_V=0.3) for i in range(5)], "id", replicates=100
    )
    assert same["intervals"]["MAE_V"] == pytest.approx([0.3, 0.3])


@pytest.mark.parametrize("value", [True, float("nan"), float("inf"), "nan"])
def test_nonfinite_and_boolean_residuals_are_rejected(value):
    with pytest.raises(ValueError):
        number(value)


def test_empty_metrics_are_blank_and_bias_is_preserved():
    assert metrics([])["MAE_V"] is None
    assert metrics([0.1, 0.3])["RMSE_V"] == pytest.approx(math.sqrt(0.05))
    assert metrics([0.1, 0.3])["bias_V"] == pytest.approx(0.2)
    assert key(dict(id="a", event="x", doi="d", source="s", where="w")) == (
        "a",
        "x",
        "d",
        "s",
        "w",
    )
