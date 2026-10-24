#!/usr/bin/env python
"""Replay fine-step SCF pairs against the predeclared diagnostic criteria.

This creates a diagnostic review only; it cannot publish an electronic or
thermal state and does not replace the preserved failed coarse-step check.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import time


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def verify(mapping):
    for path, expected in mapping.items():
        if sha(path) != expected:
            raise ValueError("Diagnostic backing source changed: " + path)


def compare(plan, records, differences):
    """Compute errors from independent energy pairs, not reported FD values."""
    acceptance = plan["acceptance"]
    reference = records["baseline"]
    coordinates, force = reference["coordinates_Bohr"], reference["gradient_Ha_Bohr"]
    checks, derivatives = [], {}
    for row in differences:
        atom, axis, step = row["atom"], row["axis"], row["step_Bohr"]
        key = (atom, axis, step)
        if key in derivatives:
            raise ValueError("Repeated displacement cannot count as independent evidence")
        pair = []
        for sign in (-1, 1):
            point = records[f"atom{atom}_axis{axis}_h{step}_sign{sign}"]
            if point["SCF_converged"] is not True or not math.isfinite(point["E_Ha"]):
                raise ValueError("Displaced SCF is unconverged or nonfinite")
            if len(point["coordinates_Bohr"]) != len(coordinates):
                raise ValueError("Displaced atom count differs")
            for i, (now, before) in enumerate(zip(point["coordinates_Bohr"], coordinates)):
                if len(now) != 3 or len(before) != 3:
                    raise ValueError("Three Cartesian coordinates required")
                for a, (v, original) in enumerate(zip(now, before)):
                    expected = original + (sign*step if (i, a) == (atom, axis) else 0.)
                    if not math.isfinite(v) or abs(v-expected) > 1e-11:
                        raise ValueError("Displaced geometry differs from the specified component")
            pair.append(point["E_Ha"])
        numeric = math.fsum([pair[1], -pair[0]])/(2*step)
        analytic = force[atom][axis]
        if abs(numeric-row["numerical"]) > 1e-12 or abs(analytic-row["analytic"]) > 1e-12:
            raise ValueError("Stored force/FD arithmetic differs")
        derivatives[key] = numeric
        checks.append(dict(atom=atom, axis=axis, step_Bohr=step, numeric=numeric,
            analytic=analytic, error_Ha_Bohr=abs(numeric-analytic)))
    pairs = sorted({(atom, axis) for atom, axis, _ in derivatives})
    fine = acceptance["finest_steps_Bohr"]
    changes = [dict(atom=a, axis=x, difference_Ha_Bohr=abs(derivatives[a, x, fine[0]]-derivatives[a, x, fine[1]]))
               for a, x in pairs]
    maximum = max(r["error_Ha_Bohr"] for r in checks)
    largest_change = max(r["difference_Ha_Bohr"] for r in changes)
    return dict(checks=checks, finest_pair_checks=changes,
        maximum_total_derivative_error_Ha_Bohr=maximum,
        maximum_finest_pair_difference_Ha_Bohr=largest_change,
        derivative_criteria_passed=(maximum <= acceptance["maximum_total_derivative_error_Ha_Bohr"]
            and largest_change <= acceptance["finest_pair_derivative_difference_Ha_Bohr"]))


def review(campaign):
    plan_path, status_path = campaign / "plan.json", campaign / "status.json"
    plan, status = (json.loads(p.read_text()) for p in (plan_path, status_path))
    if status["status"] != "fine_step_diagnostic_completed_independent_review_pending":
        raise ValueError("The live diagnostic is not complete")
    verify(plan["producer_sha256"])
    verify(plan["input_sha256"])
    verify(status["output_sha256"])
    raw = Path(plan["raw_root"]) / "diagnostic"
    evidence = json.loads((raw / "status.json").read_text())
    if evidence["status"] != "diagnostic_complete_no_promotion":
        raise ValueError("Underlying SCF diagnostic did not finish")
    settings = evidence["settings"]
    if (settings["charge"], settings["mult"], settings["grid_level"], settings["basis"], settings["xc"], settings["solvent"]) != (0, 1, 3, "def2-svp", "r2scan", "acetonitrile"):
        raise ValueError("Scientific settings differ from the declared control")
    if settings["steps_Bohr"] != [0.0001, 0.00005, 0.00001] or evidence["diagnostic_controls"]["baseline_repeats"] != 3:
        raise ValueError("Predeclared step or repetition coverage differs")
    records = {}
    for row in evidence["evaluations"]:
        p = raw / (row["name"]+".json")
        if sha(p) != row["file_sha256"] or row["name"] in records:
            raise ValueError("SCF evaluation changed or was counted twice")
        records[row["name"]] = json.loads(p.read_text())
    baseline = records["baseline"]
    flat = [(abs(v), i, a) for i, xyz in enumerate(baseline["gradient_Ha_Bohr"]) for a, v in enumerate(xyz)]
    probes = {(i, a) for _, i, a in sorted(flat, reverse=True)[:2]}
    probes.update(tuple(p) for p in evidence["requested_probes"])
    expected = {(a, x, h) for a, x in probes for h in settings["steps_Bohr"]}
    observed = {(r["atom"], r["axis"], r["step_Bohr"]) for r in evidence["differences"]}
    if expected != observed or len(records) != 3+2*len(expected):
        raise ValueError("Not all predeclared components/steps/repeats are present")
    repeats = [records[n] for n in ("baseline", "baseline_repeat1", "baseline_repeat2")]
    for r in repeats:
        if not r["SCF_converged"] or r["coordinates_Bohr"] != baseline["coordinates_Bohr"]:
            raise ValueError("Baseline repeats do not share converged geometry")
    spread = max(r["E_Ha"] for r in repeats)-min(r["E_Ha"] for r in repeats)
    gradient_spread = max(abs(v-base) for r in repeats for xyz, ref in zip(r["gradient_Ha_Bohr"], baseline["gradient_Ha_Bohr"]) for v, base in zip(xyz, ref))
    result = compare(plan, records, evidence["differences"])
    a = plan["acceptance"]
    passed = (result["derivative_criteria_passed"] and spread <= a["maximum_baseline_energy_spread_Ha"]
              and gradient_spread <= a["maximum_repeated_gradient_difference_Ha_Bohr"]
              and 0 < status["observed_peak_GPU_MiB"] <= plan["maximum_peak_GPU_MiB"])
    return dict(status="fine_step_diagnostic_independently_reviewed", checked=time.time(),
        diagnostic_criteria_passed=passed, baseline_energy_spread_Ha=spread,
        maximum_repeated_gradient_difference_Ha_Bohr=gradient_spread, **result,
        observed_peak_GPU_MiB=status["observed_peak_GPU_MiB"], SCF_evaluations=len(records),
        production_admission_changed=False, original_coarse_check_replaced=False,
        source_sha256={str(p): sha(p) for p in (plan_path, status_path, Path(__file__).resolve())},
        scope="Fine-step derivative convergence only; original state, native Hessian, G and candidate decision remain unaccepted")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("campaign", type=Path)
    p.add_argument("--out", type=Path, required=True)
    args = p.parse_args()
    result = review(args.campaign.resolve())
    with args.out.open("x") as output:
        output.write(json.dumps(result, indent=2)+"\n")
    print(json.dumps(result))


if __name__ == "__main__":
    main()
