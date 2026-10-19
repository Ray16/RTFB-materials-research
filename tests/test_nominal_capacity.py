"""Structure bookkeeping is available before quantum calculations finish."""
import csv
from types import SimpleNamespace

from redox.properties import capacity_and_proxies as capacity


def test_uncomputed_candidate_keeps_nominal_mass_and_capacity_with_no_solvation_value(tmp_path, monkeypatch):
    monkeypatch.setattr(capacity, "RESULTS", tmp_path)
    monkeypatch.setattr(capacity, "_manifest_smiles", lambda: {"uncomputed": "CC"})
    monkeypatch.setattr(capacity, "read_manifest", lambda: [dict(id="uncomputed", state="neu", charge=0),
                                                            dict(id="uncomputed", state="red1", charge=-1),
                                                            dict(id="uncomputed", state="red2", charge=-2)])
    monkeypatch.setattr(capacity, "_states", lambda gid: [])
    monkeypatch.setattr(capacity, "_sascorer", lambda: SimpleNamespace(calculateScore=lambda mol: 1.))
    capacity.main()
    rows = list(csv.DictReader((tmp_path / "capacity_and_proxies.csv").open()))
    assert len(rows) == 1 and rows[0]["n_electrons"] == "2"
    assert float(rows[0]["MW"]) > 0 and float(rows[0]["specific_capacity_mAh_g"]) > 0
    assert rows[0]["dGsolv_neutral_eV"] == ""
