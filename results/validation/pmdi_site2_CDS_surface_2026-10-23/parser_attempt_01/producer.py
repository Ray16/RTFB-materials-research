#!/usr/bin/env python
"""Isolate legacy SMD CDS step dependence without an SCF or a GPU context.

Run the installed GPU4PySCF wrapper against its CPU-only Fortran CDS library,
and separately compare the PySCF library on every sampled coordinate. This is
a diagnostic; neither the original force gate nor state admission is changed.
"""
import argparse
import ast
import ctypes
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import time
from types import SimpleNamespace

import numpy as np

from redox.core.io import atomic_write_json

STEPS = (0.001, 0.0005, 0.0001, 0.00005, 0.00001, 0.000005, 0.000001)


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def checked_displacement(reference, point, atom, axis, step, sign):
    expected = np.array(reference, dtype=float, copy=True)
    expected[atom, axis] += sign * step
    if not point.get("SCF_converged") or not np.allclose(
            expected, point["coordinates_Bohr"], atol=1e-11, rtol=0):
        raise ValueError("Stored displaced point does not match its named coordinate")


def component_diagnostic(check, cds_derivatives):
    fine = cds_derivatives[str(STEPS[-1])]
    coarse = cds_derivatives[str(check["step_Bohr"]) ]
    total_change = check["numerical"] - check["analytic"]
    cds_change = coarse - fine
    return dict(atom=check["atom"], axis=check["axis"], step_Bohr=check["step_Bohr"],
        total_numeric_minus_analytic_Ha_Bohr=total_change,
        CDS_coarse_minus_fine_Ha_Bohr=cds_change,
        remaining_difference_Ha_Bohr=total_change-cds_change)


def legacy_cpu_function(source, library):
    """Load only the pure ctypes wrapper; do not import GPU4PySCF or CuPy."""
    tree = ast.parse(source.read_text())
    function = next(n for n in tree.body
                    if isinstance(n, ast.FunctionDef) and n.name == "get_cds_legacy")
    constants = {}
    for node in tree.body:
        if isinstance(node, ast.Assign) and len(node.targets) == 1:
            target = node.targets[0]
            if isinstance(target, ast.Name) and target.id in ("solvent_db", "hartree2kcal"):
                constants[target.id] = ast.literal_eval(node.value)
    if set(constants) != {"solvent_db", "hartree2kcal"}:
        raise ValueError("Installed CDS constants changed")
    # The exact linked library is CPU-only; this interface has no CUDA calls.
    namespace = dict(constants, np=np, ctypes=ctypes, libsolvent=ctypes.CDLL(str(library)))
    exec(compile(ast.Module(body=[function], type_ignores=[]), str(source), "exec"), namespace)
    return namespace["get_cds_legacy"], constants


class CoordinateMolecule:
    """The installed CDS wrapper needs only nuclear charges and Bohr coordinates."""

    def __init__(self, charges, coordinates):
        self.natm = len(charges)
        self.charges = np.asarray(charges, dtype=int)
        self.coordinates = np.asarray(coordinates, dtype=float)

    def atom_charges(self):
        return self.charges

    def atom_coords(self, unit):
        if unit != "B":
            raise ValueError("CDS wrapper requested an unsupported coordinate unit")
        return self.coordinates


def diagnose(state, out):
    from ase.io import read
    import pyscf

    check_dir = state / "total_energy_force_check"
    summary_path = check_dir / "summary.json"
    base_path = check_dir / "baseline.json"
    geometry = state / "opt.xyz"
    summary = json.loads(summary_path.read_text())
    baseline = json.loads(base_path.read_text())
    if summary["geometry_sha256"] != digest(geometry) or not baseline["SCF_converged"]:
        raise ValueError("Stored force baseline/geometry identity is invalid")
    coordinates = np.asarray(baseline["coordinates_Bohr"], dtype=float)
    atoms = read(geometry)
    if coordinates.shape != (len(atoms), 3) or not np.allclose(
            coordinates, atoms.positions / 0.52917721067, atol=1e-11, rtol=0):
        raise ValueError("Force baseline does not use the exact stored geometry")
    site_packages = Path(pyscf.__file__).resolve().parent.parent
    gpu_root, cpu_root = site_packages / "gpu4pyscf", site_packages / "pyscf"
    sources = [gpu_root / "solvent/smd.py", cpu_root / "solvent/smd.py"]
    libraries = [gpu_root / "lib/libsolvent.so", cpu_root / "lib/libsolvent.so"]
    functions = [legacy_cpu_function(s, p) for s, p in zip(sources, libraries)]
    if functions[0][1] != functions[1][1]:
        raise ValueError("CPU/GPU wrapper solvent descriptors or conversion differ")
    molecule = CoordinateMolecule(atoms.numbers, coordinates)
    solvent = SimpleNamespace(mol=molecule, solvent="acetonitrile", solvent_descriptors=None)
    provenance = {str(p): digest(p) for p in [Path(__file__), geometry, summary_path, base_path, *sources, *libraries]}
    maximum_library_difference = 0.
    sampled_energies = []

    def energy(position, label):
        nonlocal maximum_library_difference
        molecule.coordinates = np.asarray(position, dtype=float)
        values = [float(function(solvent)[0]) for function, _ in functions]
        if not np.isfinite(values).all():
            raise ValueError("Nonfinite CDS library output")
        maximum_library_difference = max(maximum_library_difference, abs(values[0]-values[1]))
        sampled_energies.append(dict(label=label, GPU_package_CPU_library_E_Ha=values[0],
            CPU_package_library_E_Ha=values[1], coordinates_Bohr=position.tolist()))
        return values[0]

    energy(coordinates, "baseline")
    pairs = sorted({(r["atom"], r["axis"]) for r in summary["checks"]})
    records = []
    for atom, axis in pairs:
        derivatives = {}
        for step in STEPS:
            values = []
            for sign in (-1, 1):
                shifted = coordinates.copy()
                shifted[atom, axis] += sign*step
                name = f"atom{atom}_axis{axis}_h{step}_sign{sign}"
                p = check_dir / (name+".json")
                if step in (0.0005, 0.001):
                    stored = json.loads(p.read_text())
                    checked_displacement(coordinates, stored, atom, axis, step, sign)
                    provenance[str(p)] = digest(p)
                values.append(energy(shifted, name))
            derivatives[str(step)] = (values[1]-values[0])/(2*step)
        checks = [r for r in summary["checks"] if (r["atom"], r["axis"]) == (atom, axis)]
        for check in checks:
            minus = json.loads((check_dir / f"atom{atom}_axis{axis}_h{check['step_Bohr']}_sign-1.json").read_text())
            plus = json.loads((check_dir / f"atom{atom}_axis{axis}_h{check['step_Bohr']}_sign1.json").read_text())
            independent = (plus["E_Ha"]-minus["E_Ha"])/(2*check["step_Bohr"])
            if abs(independent-check["numerical"]) > 1e-12:
                raise ValueError("Stored total-energy derivative arithmetic differs")
        records.append(dict(atom=atom, axis=axis, CDS_derivative_by_step_Ha_Bohr=derivatives,
            fine_pair_difference_Ha_Bohr=abs(derivatives[str(STEPS[-1])]-derivatives[str(STEPS[-2])]),
            total_energy_comparisons=[component_diagnostic(c, derivatives) for c in checks]))
    if any(digest(p) != h for p, h in provenance.items()):
        raise ValueError("Diagnostic inputs changed during replay")
    result = dict(status="CDS_component_step_diagnostic_complete", checked=time.time(),
        atoms=len(atoms), solvent="acetonitrile", CPU_affinity=sorted(os.sched_getaffinity(0)),
        steps_Bohr=list(STEPS), SCF_calculations_performed=0, GPU_context_created=False,
        maximum_library_energy_difference_Ha=maximum_library_difference, records=records,
        source_sha256=provenance, state_admission_changed=False, original_force_gate_changed=False,
        limitations=["Only the named CDS term and sampled coordinates are isolated",
            "Fine-step CDS agreement is not full-force or Hessian validation",
            "No stationary-state energy, G, ensemble convergence or candidate admission is inferred"])
    atomic_write_json(out / "CDS_energies.json", sampled_energies)
    result["energy_samples_sha256"] = digest(out / "CDS_energies.json")
    atomic_write_json(out / "summary.json", result)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--state-dir", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    from redox.workflows.resources import allocation_guard, cpu_allocation, host_capacity

    with allocation_guard():
        host_capacity(1)
        os.sched_setaffinity(0, cpu_allocation(1))
    args.out.mkdir(parents=True, exist_ok=False)
    result = diagnose(args.state_dir.resolve(), args.out.resolve())
    print(json.dumps({k: result[k] for k in ("status", "atoms", "SCF_calculations_performed", "maximum_library_energy_difference_Ha", "records")}))


if __name__ == "__main__":
    main()
