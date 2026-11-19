"""Failures in radius and spin inputs must surface instead of switching definitions."""
import pytest


def test_sasa_failure_returns_none_not_volume_radius(monkeypatch):
    pytest.importorskip("rdkit")
    from rdkit.Chem import rdFreeSASA
    from redox.properties import solvated_reorg

    def boom(*args, **kwargs):
        raise RuntimeError("synthetic SASA failure")

    monkeypatch.setattr(rdFreeSASA, "CalcSASA", boom)
    with pytest.warns(UserWarning, match="returning None"):
        assert solvated_reorg.hydrodynamic_radius_A("c1ccccc1", probe_radius_A=1.3) is None


def test_corrupt_uma_spin_record_raises(monkeypatch, tmp_path):
    from redox.backends import dft

    gid, state = "test_only_no_such_candidate", "red1"
    record = tmp_path / gid / state / "result.json"
    record.parent.mkdir(parents=True)
    record.write_text("{not json")
    monkeypatch.setattr(dft, "UMA", tmp_path)
    with pytest.raises(RuntimeError, match="unreadable UMA spin record"):
        dft._uma_mult(gid, state, fallback=2)


def test_missing_uma_spin_record_uses_fallback(monkeypatch, tmp_path):
    from redox.backends import dft

    monkeypatch.setattr(dft, "UMA", tmp_path)
    assert dft._uma_mult("test_only_no_such_candidate", "red1", fallback=2) == 2
