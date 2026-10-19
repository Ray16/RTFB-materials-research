"""The benchmark must use the same admitted thermal model as its live Fc pair."""
import csv
import importlib.util
import json

import pytest

from redox.core.protocol import ACTIVE_SP
from redox.core.common import ROOT
from redox.validation import orop


def state(energy, thermal=.2, definition="native_solution_dft_qrrho_spin"):
    return dict(e_smd_eV=energy, g_thermal_eV=thermal, n_imag=0,
                thermal_definition=definition, thermal_qc="ok", geometry_spin_qc="ok")


def test_matching_reference_cancels_uniform_thermal_shift():
    records = [state(value) for value in (-10., -12., -50., -54.)]
    value, qc, _ = orop.matched_potential(*records)
    assert value == pytest.approx(-2.) and qc == "ok"
    for record in records:
        record["g_thermal_eV"] += .125
    assert orop.matched_potential(*records)[:2] == (value, "ok")


@pytest.mark.parametrize("change", ["legacy_molecule", "mixed_fc", "imaginary", "missing_energy", "unknown_definition"])
def test_invalid_or_mismatched_pair_has_no_referenced_prediction(change):
    records = [state(value) for value in (-10., -12., -50., -54.)]
    if change == "legacy_molecule":
        for record in records[:2]: record["thermal_definition"] = "legacy_gfn2_gas_rrho"
    elif change == "mixed_fc": records[-1]["thermal_definition"] = "legacy_gfn2_gas_rrho"
    elif change == "imaginary": records[0]["n_imag"] = 1
    elif change == "missing_energy": records[1]["e_smd_eV"] = None
    else:
        for record in records: record["thermal_definition"] = "undeclared_thermal_model"
    assert orop.matched_potential(*records)[0] is None


def test_pending_aggregation_replaces_stale_numbers_and_preserves_the_old_source(tmp_path, monkeypatch):
    monkeypatch.setattr(orop, "load_index", lambda: {1: dict(charge_ox=0, charge_red=-1, exp=-1., imp_dft=-1.1)})
    monkeypatch.setattr(orop, "state_result", lambda directory, charge:
                        (state(-10. if charge == 0 else -12., definition="legacy_gfn2_gas_rrho"), "ok"))
    monkeypatch.setattr(orop, "read_result", lambda gid, name: state(-50. if name == "ox" else -54.))
    destination = tmp_path / "orop_benchmark.csv"
    old = "sys,calc,err\n1,123.,999.\n"
    destination.write_text(old)
    rows = orop.aggregate([1], destination)
    assert rows[0]["status"] == "thermal_definition_mismatch"
    assert rows[0]["calc"] is rows[0]["err"] is None
    assert rows[0]["calibration_eligible"] is False
    saved = list(csv.DictReader(destination.open()))
    assert saved[0]["calc"] == saved[0]["err"] == ""
    history = list((tmp_path / "historical_inputs/orop_benchmark").glob("*.csv"))
    assert len(history) == 1 and history[0].read_text() == old
    orop.aggregate([], destination)
    assert list(csv.DictReader(destination.open())) == []


@pytest.mark.parametrize("change", [{}, {"reference_abs_V": 3.}, {"status": "INCOMPLETE"},
                                   {"thermal_definition": "legacy_gfn2_gas_rrho"}, {"err": "nan"}])
def test_calibration_rejects_unlabelled_historical_and_changed_reference_rows(change):
    valid = dict(status="ok", calibration_eligible="True", energy_protocol=json.dumps(ACTIVE_SP),
                 reference_abs_V=4., reference_thermal_definition="native_solution_dft_qrrho_spin",
                 thermal_definition="native_solution_dft_qrrho_spin", calc=-1., err=.1)
    assert orop.calibration_row_ok(valid, 4., "native_solution_dft_qrrho_spin")
    changed = dict(valid, **change) if change else dict(calc=-1., err=.1, charge_ox=1)
    assert not orop.calibration_row_ok(changed, 4., "native_solution_dft_qrrho_spin")


def test_unaddressed_backfilled_thermal_is_never_a_native_fallback(tmp_path, monkeypatch):
    from redox.qm import dft
    raw = dict(charge=0, mult=1, optimized=True, g_thermal_eV=None)
    (tmp_path / "result.json").write_text(json.dumps(raw))
    (tmp_path / "thermal.json").write_text(json.dumps(dict(g_thermal_eV=9.0)))
    monkeypatch.setattr(dft, "_opt_protocol_mismatch", lambda *a: "")
    monkeypatch.setattr(orop, "load_record", lambda *a: dict(e_smd_eV=-10.))
    record, _ = orop.state_result(tmp_path, 0)
    assert record["g_thermal_eV"] is None
    assert orop.matched_potential(record, state(-12.), state(-50.), state(-54.))[0] is None


def test_independent_benchmark_audit_detects_a_numeric_value_from_mismatched_models(tmp_path, monkeypatch):
    spec = importlib.util.spec_from_file_location('independent_orop_test', ROOT / 'scripts/validation/audit/recompute_axes.py')
    audit = importlib.util.module_from_spec(spec); spec.loader.exec_module(audit)
    monkeypatch.setattr(audit, 'ROOT', tmp_path)
    monkeypatch.setattr(audit, 'RES', tmp_path / 'results')
    dataset = tmp_path / 'data/raw/validation/OROP'; (dataset / 'features').mkdir(parents=True)
    (dataset / 'implicit_solvation_results.csv').write_text('exp. redox pot. [V],implicit solvation redox pot. [V]\n-2.,-2.1\n')
    (dataset / 'features/feature_input_implicit_OROP_B3LYP-D3.csv').write_text('# header\n# units\n0,0,1,0\n')
    for name, charge in [('ox',0),('red',-1)]:
        directory = tmp_path / 'calcs/orop/1' / name; directory.mkdir(parents=True)
        raw = dict(charge=charge, mult=1 if charge==0 else 2, optimized=True,
                   opt_xc='r2scan',opt_basis='def2-svp' if charge==0 else 'def2-svpd',
                   opt_disp='d4',solvent='acetonitrile')
        (directory/'result.json').write_text(json.dumps(raw))
    monkeypatch.setattr(audit,'load_record',lambda directory,*a:dict(e_smd_eV=-10. if directory.name=='ox' else -12.))
    monkeypatch.setattr(audit,'thermal_inputs',lambda *a:(0.,0,'legacy_gfn2_gas_rrho'))
    monkeypatch.setattr(audit,'state',lambda *a:dict(thermal_definition='native_solution_dft_qrrho_spin'))
    monkeypatch.setattr(audit,'fc_abs_live',lambda:4.)
    monkeypatch.setattr(audit,'read_csv',lambda *a:[dict(sys=1,charge_ox=0,exp=-2.,calc=-2.,err=0.,imp_dft=-2.1,err_impdft=-.1,
                        status='ok',thermal_definition='legacy_gfn2_gas_rrho',reference_thermal_definition='native_solution_dft_qrrho_spin',
                        reference_abs_V=4.,energy_protocol=json.dumps(ACTIVE_SP),calibration_eligible='True')])
    audit.ROWS.clear(); audit.audit_orop()
    bad={row['key'] for row in audit.ROWS if row['verdict']!='ok'}
    assert {'calc','err','status','calibration_eligible'} <= bad
