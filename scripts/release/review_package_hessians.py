#!/usr/bin/env python3
"""Reconstruct packaged mode projections and nuclear sums on a CPU.

Requires NumPy and ASE to load the separately implemented audit equations.
No xTB library, CUDA context, SCF calculation or original path is used.
This numerical replay does not repeat the original stationarity/spin/FD checks.
"""
import argparse
import importlib.util
import json
from pathlib import Path

import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("package", type=Path)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    root = args.package.resolve()
    spec = importlib.util.spec_from_file_location("package_replay", root / "code/scripts/release/reproduce_review_package.py")
    reader = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(reader)
    reader.replay(root)
    spec = importlib.util.spec_from_file_location("independent_equations", root / "code/scripts/validation/audit/native_partition.py")
    equations = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(equations)
    manifest = json.loads((root / "data/manifest.json").read_text())
    checked, geometries, frequencies = [], {}, {}
    for row in manifest["states"]:
        if row["partition_kind"] == "monatomic":
            continue
        files = row["files"]
        part = json.loads(reader.local(root / "data", files["partition"]).read_text())
        lines = reader.local(root / "data", files["geometry"]).read_text().splitlines()[2:]
        xyz = np.array([[float(v) for v in line.split()[1:]] for line in lines])
        mass = np.array(part["atomic_masses_amu"])
        if mass.shape != (len(xyz),) or not np.isfinite(mass).all() or np.any(mass <= 0):
            raise ValueError("Invalid archived masses")
        centered = xyz - np.average(xyz, axis=0, weights=mass)
        vectors = [(np.sqrt(mass)[:, None] * axis).ravel() for axis in np.eye(3)]
        vectors += [(np.sqrt(mass)[:, None] * np.cross(centered, axis)).ravel() for axis in np.eye(3)]
        rigid = []
        for vector in vectors:
            residual = vector.copy()
            for _ in range(2):
                for unit in rigid:
                    residual -= unit * np.dot(unit, residual)
            if np.linalg.norm(residual) > 1e-10 * max(1., np.linalg.norm(vector)):
                rigid.append(residual / np.linalg.norm(residual))
        basis = np.array(rigid).T
        values, vectors = np.linalg.eigh(np.eye(3*len(xyz)) - basis @ basis.T)
        internal = vectors[:, values > .5]
        hessian = np.load(reader.local(root / "data", files["hessian"]), allow_pickle=False)
        if hessian.shape != (3*len(xyz), 3*len(xyz)) or not np.isfinite(hessian).all() or abs(hessian-hessian.T).max() > 1e-6:
            raise ValueError("Invalid archived Cartesian Hessian")
        roots = np.repeat(np.sqrt(mass*(1.660539040e-27/9.10938356e-31)), 3)
        weighted = (hessian+hessian.T)/(2*np.outer(roots, roots))
        force = np.linalg.eigvalsh(internal.T @ weighted @ internal)
        modes = np.sign(force)*np.sqrt(abs(force))*219474.63067
        error = float(np.max(abs(modes-np.array(part["frequencies_cm"])), initial=0))
        if len(rigid) != part["rigid_rank"] or len(modes) != part["internal_mode_count"] or error > 1e-5:
            raise ValueError("Archived mode projection differs")
        thermal_errors = {}
        for cutoff in (25, 50, 100):
            value = equations.nuclear_equations(modes, xyz, mass, part["symmetry_number"], cutoff)
            difference = abs(value-part["nuclear_by_cutoff_eV"][str(cutoff)])
            if difference > 1e-8:
                raise ValueError("Archived nuclear partition differs")
            thermal_errors[str(cutoff)] = difference
        checked.append(dict(id=row["id"], frequency_max_error_cm=error, nuclear_errors_eV=thermal_errors))
        geometries[row["id"]] = xyz
        frequencies[row["id"]] = modes
    restart_checks = []
    if "auxiliary" in manifest:
        package = json.loads((root / "package_manifest.json").read_text())
        epoch = reader.local(root, package["evidence_epoch"])
        controls = json.loads((epoch / "quantum_restart.json").read_text())
        for control in controls["controls"]:
            a, b = [manifest["aliases"]["restart:" + control["case"] + "__" + phase] for phase in ("baseline", "resumed")]
            x, y = geometries[a], geometries[b]
            x, y = x-x.mean(axis=0), y-y.mean(axis=0)
            u, _, vt = np.linalg.svd(x.T @ y)
            rotation = u @ np.diag([1, 1, np.linalg.det(u @ vt)]) @ vt
            rmsd = float(np.sqrt(np.mean(np.sum((x @ rotation-y)**2, axis=1))))
            difference = float(np.max(abs(np.sort(frequencies[a])-np.sort(frequencies[b]))))
            if (abs(rmsd-control["differences"]["aligned_RMSD_Angstrom"]) > 1e-8
                    or abs(difference-control["differences"]["frequency_cm"]) > 1e-8):
                raise ValueError("Archived restart geometry/frequency comparison differs")
            restart_checks.append(dict(case=control["case"], aligned_RMSD_Angstrom=rmsd, frequency_cm=difference))
    result = dict(status="packaged_Hessians_modes_and_nuclear_sums_reconstructed", checked=len(checked), states=checked,
        restart_geometry_frequency_checks=restart_checks,
        original_cluster_paths_opened=False, quantum_calculations_run=False, xTB_library_loaded=False,
        original_physical_QC_repeated=False, experimental_accuracy_inferred=False)
    if args.out:
        with args.out.open("x") as stream:
            json.dump(result, stream, indent=2)
            stream.write("\n")
    print(json.dumps({k:v for k,v in result.items() if k != "states"}))


if __name__ == "__main__":
    main()
