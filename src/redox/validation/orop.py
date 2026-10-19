"""OROP benchmark with live, explicitly matched molecular and Fc thermochemistry.

Experimental and original-dataset predictions remain source measurements. Current
predictions require admitted active SPs and the same selected thermal definition
for both molecular states and both Fc states. Pending rows are always written.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import shutil
import time

from redox.core.common import ROOT, RESULTS, DFT, apply_thermal_policy, free_energy, read_result
from redox.core.io import atomic_text_writer, atomic_write_json
from redox.core.protocol import ACTIVE_SP, load_record, record_path
from redox.properties.potentials import thermal_definition

OROP = ROOT / "data/raw/validation/OROP"
XYZ = OROP / "implicit_optimized_xyz"
CALC = ROOT / "calcs/orop"
DEFAULT = [103, 96, 101, 90, 16, 115, 104, 105, 109, 108, 171, 107, 192, 191]
COLUMNS = ["sys", "charge_ox", "exp", "calc", "imp_dft", "err", "err_impdft",
           "status", "thermal_definition", "reference_thermal_definition",
           "reference_abs_V", "energy_protocol", "calibration_eligible", "qc_detail"]


def load_index():
    import numpy as np
    import pandas as pd
    features = np.loadtxt(OROP / "features/feature_input_implicit_OROP_B3LYP-D3.csv",
                          skiprows=2, delimiter=",", ndmin=2)
    reference = pd.read_csv(OROP / "implicit_solvation_results.csv")
    reference.columns = [column.strip() for column in reference.columns]
    index = {}
    for offset, feature in enumerate(features):
        system = offset + 1
        ox, red = XYZ / f"{system}-1-MeCN.xyz", XYZ / f"{system}-2-MeCN.xyz"
        if not (ox.is_file() and red.is_file()):
            continue
        charge, mult = int(feature[3]), int(feature[2])
        index[system] = dict(charge_ox=charge, mult_ox=mult, charge_red=charge - 1,
                             mult_red=1 if mult == 2 else 2, xyz_ox=ox, xyz_red=red,
                             exp=float(reference["exp. redox pot. [V]"].iloc[offset]),
                             imp_dft=float(reference["implicit solvation redox pot. [V]"].iloc[offset]))
    return index


def matched_potential(ox, red, fc_ox, fc_neu):
    """Pure admission/arithmetic, with a reason for each unavailable prediction."""
    records = (ox, red, fc_ox, fc_neu)
    values = [free_energy(record) for record in records]
    names = ("ox", "red", "Fc_ox", "Fc_neu")
    bad = [f"{name}:{qc if qc != 'ok' else 'missing_energy'}"
           for name, (energy, qc) in zip(names, values) if energy is None or qc != "ok"]
    if bad:
        return None, "INCOMPLETE", ";".join(bad)
    definitions = {thermal_definition(record) for record in records}
    if not definitions <= {"legacy_gfn2_gas_rrho", "native_solution_dft_qrrho_spin"}:
        return None, "invalid_thermal_definition", repr(sorted(definitions, key=str))
    if len(definitions) != 1:
        return None, "thermal_definition_mismatch", ";".join(
            f"{name}:{thermal_definition(record)}" for name, record in zip(names, records))
    go, gr, gfo, gfn = [value for value, _ in values]
    return -(gr - go) - (gfo - gfn), "ok", "matched admitted free energies"


def calibration_row_ok(row, reference_abs, reference_definition):
    """Unlabelled historical or stale-reference rows cannot calibrate the live model."""
    from redox.core.common import to_float
    try:
        protocol = json.loads(row.get("energy_protocol", ""))
    except (ValueError, TypeError):
        return False
    reference = to_float(row.get("reference_abs_V"))
    return (row.get("status") == "ok"
            and str(row.get("calibration_eligible")).lower() == "true"
            and protocol == ACTIVE_SP
            and reference_abs is not None and reference is not None
            and abs(reference - reference_abs) < 1e-8
            and row.get("thermal_definition") == row.get("reference_thermal_definition") == reference_definition
            and to_float(row.get("calc")) is not None and to_float(row.get("err")) is not None)


def state_result(directory, expected_charge):
    """Read one external benchmark state through the production energy/thermal gates."""
    try:
        raw = json.loads((directory / "result.json").read_text())
        if raw["charge"] != expected_charge:
            return None, "charge_mismatch"
        from redox.qm.dft import _opt_protocol_mismatch
        if _opt_protocol_mismatch(raw, expected_charge):
            return None, "stale_geometry_protocol"
        energy = load_record(directory, raw["charge"], raw["mult"])
        if energy is None:
            return None, "missing_or_invalid_active_sp"
        selected = apply_thermal_policy(dict(raw, e_smd_eV=energy["e_smd_eV"]), directory)
        return selected, "ok"
    except (OSError, ValueError, KeyError, TypeError):
        return None, "missing_or_malformed_result"


def aggregate(systems, output=None):
    index = load_index()
    fc_ox, fc_neu = read_result("ferrocene", "ox"), read_result("ferrocene", "neu")
    gfo, qfo = free_energy(fc_ox)
    gfn, qfn = free_energy(fc_neu)
    reference = gfo - gfn if (gfo is not None and gfn is not None and qfo == qfn == "ok"
                            and thermal_definition(fc_ox) == thermal_definition(fc_neu)) else None
    rows = []
    for system in systems:
        source = index.get(system)
        if source is None:
            raise ValueError(f"OROP system {system} has no source index/geometry")
        ox, qo = state_result(CALC / str(system) / "ox", source["charge_ox"])
        red, qr = state_result(CALC / str(system) / "red", source["charge_red"])
        value, status, detail = matched_potential(ox, red, fc_ox, fc_neu)
        if qo != "ok" or qr != "ok":
            detail = f"ox:{qo};red:{qr};{detail}"
        rows.append(dict(sys=system, charge_ox=source["charge_ox"], exp=source["exp"],
                         calc=round(value, 3) if value is not None else None,
                         imp_dft=round(source["imp_dft"], 3),
                         err=round(value - source["exp"], 3) if value is not None else None,
                         err_impdft=round(source["imp_dft"] - source["exp"], 3), status=status,
                         thermal_definition=thermal_definition(ox) if thermal_definition(ox) == thermal_definition(red) else "mismatch",
                         reference_thermal_definition=thermal_definition(fc_ox), reference_abs_V=reference,
                         energy_protocol=json.dumps(ACTIVE_SP, sort_keys=True),
                         calibration_eligible=status == "ok", qc_detail=detail))
    output = Path(output) if output is not None else RESULTS / "orop_benchmark.csv"
    # Preserve any earlier derived benchmark, including its former reference treatment.
    if output.is_file():
        digest = hashlib.sha256(output.read_bytes()).hexdigest()
        history = output.parent / "historical_inputs/orop_benchmark" / f"{digest}.csv"
        history.parent.mkdir(parents=True, exist_ok=True)
        if not history.exists():
            shutil.copy2(output, history)
    with atomic_text_writer(output, newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=COLUMNS)
        writer.writeheader()
        writer.writerows(rows)
    print(f"[orop] {sum(row['calibration_eligible'] for row in rows)}/{len(rows)} matched predictions; {output}")
    return rows


def aggregate_all():
    systems = sorted(int(path.name) for path in CALC.iterdir() if path.is_dir() and path.name.isdigit())
    return aggregate(systems)


def run_system(system, backend="gpu", force=False):
    """Compute geometry, active SP and native thermal records with explicit failure."""
    if backend != "gpu":
        raise ValueError("Native benchmark thermochemistry requires the validated GPU backend")
    from redox.qm import dft, sp, native_thermal
    source = load_index()[system]
    failures = 0
    for state in ("ox", "red"):
        directory = CALC / str(system) / state
        raw_path = directory / "result.json"
        charge, mult = source[f"charge_{state}"], source[f"mult_{state}"]
        try:
            if raw_path.exists() and force:
                history = ROOT / "results/history/orop" / str(system) / state / str(time.time_ns())
                shutil.copytree(directory, history)
            if not raw_path.exists() or force:
                directory.mkdir(parents=True, exist_ok=True)
                result = dft.dft_smd(source[f"xyz_{state}"], charge, mult, do_opt=True,
                                     do_freq=False, backend=backend, opt_out=directory / "opt.xyz")
                atomic_write_json(raw_path, result)
            raw = json.loads(raw_path.read_text())
            if (raw.get("charge") != charge or dft._opt_protocol_mismatch(raw, charge)
                    or raw.get("converged_smd") is not True):
                raise RuntimeError("Existing benchmark geometry requires separate repair")
            sp.ensure(f"orop{system}", state, directory, charge, raw["mult"], backend=backend)
            native_thermal.select_policy(directory)
            native_thermal.ensure_native_thermal(f"orop{system}", state, state_dir=directory)
            print(f"[orop] {system}/{state}: native records complete", flush=True)
        except Exception as exc:
            failures += 1
            print(f"[orop] {system}/{state}: pending/failed: {exc!r}", flush=True)
    return int(bool(failures))


def fanout(systems, force=False):
    raise RuntimeError("Legacy unreserved OROP fanout is retired; use redox.workflows.gpu")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--only", type=int)
    mode.add_argument("--aggregate")
    mode.add_argument("--fanout")
    parser.add_argument("--out", type=Path, help="benchmark CSV destination")
    parser.add_argument("--backend", default="gpu", choices=("gpu", "cpu"))
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args(argv)
    if args.fanout:
        fanout([])
    if args.only is not None:
        return run_system(args.only, args.backend, args.force)
    if args.aggregate == "all":
        systems = sorted(int(path.name) for path in CALC.iterdir() if path.is_dir() and path.name.isdigit())
    else:
        systems = DEFAULT if args.aggregate == "default" else [int(value) for value in args.aggregate.split(",")]
    aggregate(systems, args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
