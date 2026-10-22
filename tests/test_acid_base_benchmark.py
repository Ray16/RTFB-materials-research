import importlib.util
from pathlib import Path

import pytest


def module():
    path = Path(__file__).resolve().parents[1] / "scripts/validation/benchmark_acid_base.py"
    spec = importlib.util.spec_from_file_location("acid_base_term_check", path)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


def test_capture_sign_and_temperature_are_defined():
    code = module()
    assert code.capture_deltaG(10, 10) == 0
    assert code.capture_deltaG(10, 11) == pytest.approx(-5.708, abs=.001)
    assert code.capture_deltaG(11, 10) > 0
    with pytest.raises(ValueError):
        code.capture_deltaG(10, 11, 0)


def test_exact_graph_matching_does_not_substitute_a_related_amine():
    code = module()
    table = code.BasicityTable()
    table.feed('<table><tr><td>108-18-9</td><td>Diisopropylamine</td><td></td>'
               '<td><a href="https://example.org/original">18.82</a></td><td></td>'
               '<td>CC(C)NC(C)C</td></tr></table>')
    assert code.exact_base(table, "CCN(C(C)C)C(C)C") is None
    result = code.exact_base(table, "CC(C)NC(C)C")
    assert result["pKaH_MeCN"] == 18.82
    assert result["source_urls"] == ["https://example.org/original"]
