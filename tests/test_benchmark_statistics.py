"""Sparse or invalid benchmarks cannot fabricate numerical uncertainty."""
import importlib.util
import json
import math

import pytest

from redox.core.common import ROOT


@pytest.fixture
def scorer():
    spec = importlib.util.spec_from_file_location("benchmark_statistics_test", ROOT / "scripts/validation/benchmark/score_benchmark_mecn.py")
    module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module


@pytest.mark.parametrize("residuals", [[], [.2]])
def test_unestimable_sample_sd_is_blank_and_serializes_as_json_null(scorer, residuals):
    report = scorer.stats(residuals)
    assert report["n"] == len(residuals)
    assert report["SD"] is None
    json.dumps(report, allow_nan=False)
    if residuals:
        assert report["MAE"] == report["RMSE"] == report["bias"] == .2
    else:
        assert report["MAE"] is report["RMSE"] is report["bias"] is None


def test_observed_bias_is_retained_in_rmse(scorer):
    report = scorer.stats([.1, .3])
    assert report["bias"] == pytest.approx(.2)
    assert report["RMSE"] == pytest.approx(math.sqrt(.05))
    assert report["SD"] == pytest.approx(math.sqrt(.02))


@pytest.mark.parametrize("value", [float("nan"), float("inf"), -float("inf")])
def test_nonfinite_values_cannot_enter_benchmark_statistics(scorer, value):
    assert scorer._f(str(value)) is None
    with pytest.raises(ValueError):
        scorer.stats([value])
