#!/usr/bin/env python
"""Independently reconstruct the CDS diagnostic from its stored energy pairs."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import time


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def review(folder):
    summary_path, energies_path = folder / "summary.json", folder / "CDS_energies.json"
    summary = json.loads(summary_path.read_text())
    for path, expected in summary["source_sha256"].items():
        if sha(path) != expected:
            raise ValueError("CDS source bytes changed: " + path)
    if sha(energies_path) != summary["energy_samples_sha256"]:
        raise ValueError("CDS samples changed")
    rows = json.loads(energies_path.read_text())
    samples = {r["label"]: r for r in rows}
    if len(samples) != len(rows) or len(rows) != 1 + 2*len(summary["steps_Bohr"])*len(summary["records"]):
        raise ValueError("CDS sample inventory differs")
    reference = samples["baseline"]["coordinates_Bohr"]
    comparisons, largest_remaining = [], 0.
    library_differences = [abs(r["GPU_package_CPU_library_E_Ha"]-r["CPU_package_library_E_Ha"]) for r in rows]
    for record in summary["records"]:
        atom, axis = record["atom"], record["axis"]
        derivatives = {}
        for step in summary["steps_Bohr"]:
            pair = []
            for sign in (-1, 1):
                point = samples[f"atom{atom}_axis{axis}_h{step}_sign{sign}"]
                coords = point["coordinates_Bohr"]
                if len(coords) != len(reference):
                    raise ValueError("CDS atom count differs")
                for i, (current, base) in enumerate(zip(coords, reference)):
                    for a, (v, original) in enumerate(zip(current, base)):
                        expected = original + (sign*step if (i, a) == (atom, axis) else 0.)
                        if not math.isfinite(v) or abs(v-expected) > 1e-11:
                            raise ValueError("CDS coordinate changed outside the specified displacement")
                pair.append(point["GPU_package_CPU_library_E_Ha"])
            derivative = math.fsum([pair[1], -pair[0]])/(2*step)
            if abs(derivative-record["CDS_derivative_by_step_Ha_Bohr"][str(step)]) > 1e-13:
                raise ValueError("CDS derivative arithmetic differs")
            derivatives[step] = derivative
        fine = derivatives[min(derivatives)]
        source = next(Path(p) for p in summary["source_sha256"] if p.endswith("total_energy_force_check/summary.json"))
        original = json.loads(source.read_text())
        for step in (0.0005, 0.001):
            stored = next(r for r in original["checks"] if (r["atom"], r["axis"], r["step_Bohr"]) == (atom, axis, step))
            pair = [json.loads((source.parent/f"atom{atom}_axis{axis}_h{step}_sign{s}.json").read_text())["E_Ha"] for s in (-1, 1)]
            total_fd = math.fsum([pair[1], -pair[0]])/(2*step)
            remaining = total_fd-stored["analytic"]-derivatives[step]+fine
            reported = next(r for r in record["total_energy_comparisons"] if r["step_Bohr"] == step)
            if abs(remaining-reported["remaining_difference_Ha_Bohr"]) > 1e-13:
                raise ValueError("CDS isolation does not reconstruct")
            largest_remaining = max(largest_remaining, abs(remaining))
            comparisons.append(dict(atom=atom, axis=axis, step_Bohr=step,
                remaining_difference_Ha_Bohr=remaining))
    if not math.isclose(max(library_differences), summary["maximum_library_energy_difference_Ha"], rel_tol=0, abs_tol=1e-15):
        raise ValueError("Library comparison differs")
    return dict(status="CDS_energy_pairs_and_total_force_isolation_independently_reconstructed",
        checked=time.time(), energy_coordinates_checked=len(rows), CDS_derivatives_reconstructed=len(rows)//2,
        total_force_comparisons=comparisons, maximum_remaining_difference_Ha_Bohr=largest_remaining,
        maximum_library_energy_difference_Ha=max(library_differences),
        source_sha256={str(p): sha(p) for p in (summary_path, energies_path, Path(__file__).resolve())},
        production_admission_changed=False,
        scope="Stored energy/coordinate and component-isolation replay; no new SCF, full-force acceptance, Hessian or G")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("folder", type=Path)
    p.add_argument("--out", type=Path, required=True)
    args = p.parse_args()
    result = review(args.folder.resolve())
    with args.out.open("x") as output:
        output.write(json.dumps(result, indent=2)+"\n")
    print(json.dumps(result))


if __name__ == "__main__":
    main()
