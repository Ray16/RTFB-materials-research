#!/usr/bin/env python
"""Compare conservative SMD gradients with energy differences at a saved geometry.

Grid refinements are diagnostics; no geometry or production record is promoted.
Each SCF is checkpointed on shared storage before the next displacement.
"""

import argparse
from contextlib import nullcontext
import json
import os
from pathlib import Path
import socket
import time

import numpy as np

from redox.core.io import atomic_write_json
from redox.workflows.chunked import file_hash


def parse_probe(value):
    atom, axis = map(int, value.split(":"))
    if atom < 0 or axis not in (0, 1, 2):
        raise argparse.ArgumentTypeError("Use a nonnegative atom and axis 0, 1 or 2")
    return atom, axis


def reusable_record(parent, name, coordinates):
    path = parent / (name + ".json")
    if not path.exists():
        return None
    record = json.loads(path.read_text())
    if (
        record.get("SCF_converged") is not True
        or not np.isfinite(record.get("E_Ha", np.nan))
        or not np.array_equal(np.asarray(record["coordinates_Bohr"]), coordinates)
    ):
        raise ValueError(
            "Previous diagnostic evaluation has different coordinates or QC"
        )
    return dict(record, reused_from=str(path.resolve()), reused_sha256=file_hash(path))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--xyz", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--grid-level", type=int, choices=(3, 5), required=True)
    parser.add_argument("--probe", type=parse_probe, action="append", default=[])
    parser.add_argument("--reuse", type=Path)
    parser.add_argument("--steps", type=float, nargs="+", default=[.0005, .001])
    parser.add_argument("--baseline-repeats", type=int, default=1)
    parser.add_argument("--top-components", type=int, default=1)
    parser.add_argument("--omit-reference-axis", action="store_true")
    parser.add_argument("--conservative-cds", action="store_true",
                        help="Use the independently qualified adaptive CDS energy derivative")
    args = parser.parse_args()
    if (not all(np.isfinite(h) and h > 0 for h in args.steps)
            or len(set(args.steps)) != len(args.steps)
            or args.baseline_repeats < 1 or args.top_components < 1):
        parser.error("Distinct positive derivative steps and positive repeat/component counts required")
    from redox.workflows.resources import verify_gate_owner

    verify_gate_owner()
    if len(os.sched_getaffinity(0)) != 2:
        raise RuntimeError("Two disjoint CPU slots required")
    from ase.io import read
    from pyscf import gto
    from redox.qm.dft import _build_smd_mf, _free_gpu_pool
    from redox.qm.stable_geometry_scf import stabilize

    args.out.mkdir(parents=True, exist_ok=False)
    atoms = read(args.xyz)
    symbols = atoms.get_chemical_symbols()
    coordinates = atoms.positions / 0.52917721067
    settings = dict(
        xc="r2scan",
        disp="d4",
        basis="def2-svp",
        solvent="acetonitrile",
        solvent_model="SMD",
        grid_level=args.grid_level,
        grid_response=True,
        conv_tol=1e-12,
        conv_tol_grad=1e-8,
        charge=0,
        mult=1,
        steps_Bohr=args.steps,
    )
    state = dict(
        status="computing",
        host=socket.gethostname(),
        pid=os.getpid(),
        process_start=Path("/proc/self/stat").read_text().rsplit(")", 1)[1].split()[19],
        cpus=sorted(os.sched_getaffinity(0)),
        gpu=os.environ.get("CUDA_VISIBLE_DEVICES"),
        source_xyz=str(args.xyz.resolve()),
        source_geometry_sha256=file_hash(args.xyz),
        settings=settings,
        producer_sha256=file_hash(__file__),
        evaluations=[],
        started=time.time(),
    )
    atomic_write_json(args.out / "status.json", state)
    (args.out / "producer_source.py").write_bytes(Path(__file__).read_bytes())
    force_context = None
    if args.conservative_cds:
        from redox.core.common import ROOT
        from redox.workflows import polymer_campaign as api
        force_path = ROOT / "scripts/pipeline/resume_polymer_forces_adaptive.py"
        module = api.load_module("diagnostic_adaptive_CDS", force_path)
        force_context = module.adaptive_forces
        settings["CDS_gradient"] = "adaptive_energy_derivative_tolerance_1e-6"
        state["force_producer_sha256"] = {str(p): file_hash(p) for p in
            [force_path, ROOT / "scripts/pipeline/resume_polymer_forces.py", ROOT / "src/redox/qm/cds_finite.py"]}
    if args.top_components > 3 * len(atoms):
        raise ValueError("More force probes than Cartesian components")
    state["diagnostic_controls"] = dict(baseline_repeats=args.baseline_repeats,
        top_components=args.top_components, omit_reference_axis=args.omit_reference_axis)
    baseline = None
    reuse_ready = False
    if any(atom >= len(atoms) for atom, _ in args.probe):
        raise ValueError("Probe atom outside this geometry")
    if args.reuse:
        parent = json.loads((args.reuse / "status.json").read_text())
        if (
            parent["status"] != "diagnostic_complete_no_promotion"
            or parent["settings"] != settings
            or parent["source_geometry_sha256"] != state["source_geometry_sha256"]
        ):
            raise ValueError("Previous diagnostic protocol or geometry differs")
        state["parent_status_sha256"] = file_hash(args.reuse / "status.json")
    state["requested_probes"] = args.probe

    def solve(position, name, gradient=False):
        if reuse_ready and not gradient:
            previous = reusable_record(args.reuse, name, position)
            if previous:
                atomic_write_json(args.out / (name + ".json"), previous)
                state["evaluations"].append(
                    dict(
                        name=name,
                        file_sha256=file_hash(args.out / (name + ".json")),
                        reused=True,
                    )
                )
                atomic_write_json(args.out / "status.json", state)
                return previous["E_Ha"], None, None
        mol = gto.M(
            atom=list(zip(symbols, position.tolist())),
            unit="Bohr",
            charge=0,
            spin=0,
            basis="def2-svp",
            verbose=0,
            max_memory=6000,
        )
        mf = _build_smd_mf(mol, "r2scan", "acetonitrile", "gpu", disp="d4")
        mf.grids.level = args.grid_level
        mf.conv_tol, mf.conv_tol_grad, mf.max_cycle = 1e-12, 1e-8, 300
        stabilize(mf, args.out / "scf_recovery.jsonl")
        mf.kernel(dm0=baseline)
        if not mf.converged:
            raise RuntimeError("Diagnostic SCF did not converge")
        record = dict(
            name=name,
            E_Ha=float(mf.e_tot),
            SCF_converged=True,
            coordinates_Bohr=position.tolist(),
        )
        raw_density = mf.make_rdm1()
        density = np.asarray(
            raw_density.get() if hasattr(raw_density, "get") else raw_density
        )
        grad = None
        if gradient:
            method = mf.nuc_grad_method()
            method.grid_response = True
            context = force_context(args.out / (name + "_cds.jsonl")) if force_context else nullcontext()
            with context:
                raw_gradient = method.kernel()
            grad = np.asarray(
                raw_gradient.get() if hasattr(raw_gradient, "get") else raw_gradient
            )
            record.update(
                gradient_Ha_Bohr=grad.tolist(),
                maximum_gradient_component=float(abs(grad).max()),
            )
        atomic_write_json(args.out / (name + ".json"), record)
        state["evaluations"].append(
            dict(name=name, file_sha256=file_hash(args.out / (name + ".json")))
        )
        atomic_write_json(args.out / "status.json", state)
        del mf, mol
        _free_gpu_pool()
        return record["E_Ha"], density, grad

    try:
        energy, baseline, gradient = solve(coordinates, "baseline", True)
        repeats = []
        for repeat in range(1, args.baseline_repeats):
            again, _, repeated_gradient = solve(coordinates, "baseline_repeat" + str(repeat), True)
            repeats.append(dict(repeat=repeat, energy_difference_Ha=again-energy,
                maximum_gradient_difference_Ha_Bohr=float(np.max(abs(repeated_gradient-gradient)))))
        state["baseline_repeat_checks"] = repeats
        atomic_write_json(args.out / "status.json", state)
        if args.reuse:
            previous = reusable_record(args.reuse, "baseline", coordinates)
            if (
                abs(previous["E_Ha"] - energy) > 1e-8
                or np.max(abs(np.asarray(previous["gradient_Ha_Bohr"]) - gradient))
                > 1e-8
            ):
                raise ValueError(
                    "Fresh baseline does not replay the previous SCF/gradient"
                )
            reuse_ready = True
            state["baseline_replayed"] = True
        worst = tuple(
            int(x) for x in np.unravel_index(abs(gradient).argmax(), gradient.shape)
        )
        top = [tuple(map(int, np.unravel_index(index, gradient.shape)))
               for index in np.argsort(abs(gradient).ravel())[-args.top_components:][::-1]]
        axes = list(dict.fromkeys([*top, *([] if args.omit_reference_axis else [(0, 1)]), *args.probe]))
        differences = []
        for atom, axis in axes:
            for step in settings["steps_Bohr"]:
                energies = []
                for sign in (-1, 1):
                    shifted = coordinates.copy()
                    shifted[atom, axis] += sign * step
                    energy, _, _ = solve(
                        shifted, f"atom{atom}_axis{axis}_h{step}_sign{sign}"
                    )
                    energies.append(energy)
                numerical = (energies[1] - energies[0]) / (2 * step)
                differences.append(
                    dict(
                        atom=atom,
                        axis=axis,
                        step_Bohr=step,
                        analytic=float(gradient[atom, axis]),
                        numerical=numerical,
                        absolute_difference=abs(numerical - gradient[atom, axis]),
                    )
                )
        state.update(
            status="diagnostic_complete_no_promotion",
            worst_component=worst,
            differences=differences,
            finished=time.time(),
        )
        atomic_write_json(args.out / "status.json", state)
    except BaseException as error:
        state.update(status="failed_preserved", error=repr(error), finished=time.time())
        atomic_write_json(args.out / "status.json", state)
        raise


if __name__ == "__main__":
    main()
