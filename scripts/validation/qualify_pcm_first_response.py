#!/usr/bin/env python
"""Test blocked PCM response on an exact saved polymer geometry.

Compare the analytic MO/occupied potential derivative to fixed-density finite
differences. This is an isolated response/resource check, not a native Hessian,
free energy, optimizer minimum or candidate-admission certificate.
"""
import argparse
from contextlib import redirect_stderr, redirect_stdout
import gc
import json
import os
from pathlib import Path
import time

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api

ADAPTER = ROOT / "scripts/pipeline/smd_first_response_workspace.py"
TENSORS = ROOT / "scripts/pipeline/pcm_first_derivative_blocks.py"
THREADS = ("OMP_NUM_THREADS", "OMP_THREAD_LIMIT", "MKL_NUM_THREADS",
           "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS", "XTB_THREADS")


def checked_controls(folder):
    sources = {}
    for mult in (1, 2):
        path = folder / "controls" / str(mult) / "status.json"
        status = json.loads(path.read_text())
        if status["status"] != "full_original_compact_equivalence_passed" or status["mult"] != mult:
            raise ValueError("Both full spin controls must pass")
        api.verify_sources(status["source_sha256"])
        api.verify_sources(status["evidence_sha256"])
        stats_path = path.parent / "compact/execution/first_derivatives.json"
        stats = json.loads(stats_path.read_text())
        execution = json.loads((path.parent / "compact/execution/status.json").read_text())
        if (stats["factories"] < 1 or stats["primitive_calls"] < 1 or stats["rectangles"] < 1
                or stats["maximum_primitive_points"] > 2*stats["tile_points"]
                or stats["maximum_output_bytes"] > stats["maximum_output_bound_bytes"]
                or not stats["unchanged_pairwise_primitive"]
                or execution["budget_GiB"] != 28
                or execution.get("PCM_first_derivative_equations_changed") is not False
                or execution["source_sha256"].get(str(ADAPTER)) != api.digest(ADAPTER)):
            raise ValueError("Bounded first derivatives were not exercised under the required profile")
        api.verify_sources(stats["source_sha256"])
        sources.update({str(path): api.digest(path), str(stats_path): api.digest(stats_path)})
        sources.update(status["source_sha256"])
    return sources


def prepare(args):
    source = args.source.resolve()
    reader = api.load_module("response_checked_energy_parent", ROOT / "scripts/pipeline/workspace_polymer_campaign.py")
    record, parent_path = reader.checked_energy(source)
    if (record["task"]["charge"], record["task"]["mult"]) != (0, 1):
        raise ValueError("This response pilot is scoped to an existing neutral singlet")
    controls = checked_controls(args.controls.resolve())
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    xyz = out / "input.xyz"
    original = Path(record["state_dir"]) / "opt.xyz"
    xyz.write_bytes(original.read_bytes())
    plan = dict(task=record["task"], source_campaign=str(source),
        control_root=str(args.controls.resolve()), input_xyz=str(xyz),
        raw_root=str(ROOT / "calcs/validation" / out.name),
        input_sha256={str(p): api.digest(p) for p in (source / "plan.json", parent_path, original, xyz)},
        producer_sha256=dict(controls, **{str(p): api.digest(p) for p in
            (Path(__file__).resolve(), ADAPTER, TENSORS)}),
        settings=dict(xc="r2scan", dispersion="d4", basis="def2-svp", solvent="acetonitrile",
            grid_level=3, conv_tol=1e-12, conv_tol_grad=1e-8, max_cycle=300,
            response="fixed AO density and baseline MO coefficients; electrostatic solvent potential only",
            steps_Bohr=[0.0005, 0.001]),
        acceptance=dict(maximum_response_relative_error=.01, maximum_GPU_MiB=32768,
            minimum_GPU_MiB=31744, minimum_host_RAM_GiB=96), workspace_GiB=28,
        native_Hessian_qualified=False, G_computed=False, prepared=time.time())
    atomic_write_json(out / "plan.json", plan)
    print(json.dumps(dict(status="exact_geometry_response_pilot_prepared", atoms=record["task"]["atoms"])))


def run(args):
    from redox.workflows.resources import verify_gate_owner, host_capacity, gpu_suitability
    verify_gate_owner()
    if len(os.sched_getaffinity(0)) != 2 or any(os.environ.get(k) != "2" for k in THREADS):
        raise ValueError("Two disjoint CPU cores and all native caps=2 required")
    out = args.out.resolve()
    plan = api.read_plan(out)
    api.verify_sources(plan["producer_sha256"])
    api.verify_sources(plan["input_sha256"])
    checked_controls(Path(plan["control_root"]))
    host_capacity(2, minimum_ram_gib=plan["acceptance"]["minimum_host_RAM_GiB"])
    if gpu_suitability(int(os.environ["CUDA_VISIBLE_DEVICES"]))["memory_MiB"] < plan["acceptance"]["minimum_GPU_MiB"]:
        raise ValueError("Device is below the declared response test class")
    raw = Path(plan["raw_root"])
    raw.mkdir(parents=True, exist_ok=False)
    status_path = out / "status.json"
    if status_path.exists():
        raise FileExistsError("Preserve the earlier response attempt")
    status = dict(status="response_pilot_running", started=time.time(), host=os.uname().nodename,
        pid=os.getpid(), process_start=Path("/proc/self/stat").read_text().rsplit(")",1)[1].split()[19],
        CPU_affinity=sorted(os.sched_getaffinity(0)), GPU=os.environ["CUDA_VISIBLE_DEVICES"],
        native_Hessian_qualified=False, G_computed=False, candidate_admission=False)
    atomic_write_json(status_path, status)
    monitor = None
    try:
        from ase.io import read
        import numpy as np
        import cupy
        from pyscf import gto
        from gpu4pyscf.solvent.hessian import pcm
        from redox.qm import dft
        from redox.qm.stable_geometry_scf import stabilize
        atoms = read(plan["input_xyz"])
        coordinates = atoms.positions / 0.52917721067
        probes = [(0, 1), (next((i for i, symbol in enumerate(atoms.symbols) if symbol == "O"), len(atoms)-1), 2)]
        adapter = api.load_module("qualified_first_response_adapter", ADAPTER)
        engine = api.load_module("first_response_memory_observer", ROOT / "scripts/pipeline/run_pi_stacking.py")

        def mean_field(coords):
            mol = gto.M(atom=list(zip(atoms.get_chemical_symbols(), coords.tolist())), unit="Bohr",
                basis="def2-svp", charge=0, spin=0, verbose=0, max_memory=6000)
            mf = dft._build_smd_mf(mol, "r2scan", "acetonitrile", "gpu", disp="d4")
            mf.grids.level = 3
            mf.conv_tol, mf.conv_tol_grad, mf.max_cycle = 1e-12, 1e-8, 300
            return mf

        with engine.MemoryMonitor() as monitor, (raw / "calculation.log").open("x") as log, redirect_stdout(log), redirect_stderr(log):
            with adapter.compact_smd_execution(raw / "execution", budget_GiB=28):
                mf = mean_field(coordinates)
                stabilize(mf, raw / "SCF_recovery.jsonl")
                mf.chkfile = str(raw / "baseline_SCF.h5")
                mf.kernel()
                if not mf.converged:
                    raise ValueError("Response baseline SCF did not converge")
                dm, mo, occupied = mf.make_rdm1(), mf.mo_coeff, mf.mo_occ
                baseline = dict(E_Ha=float(mf.e_tot), coordinates_Bohr=coordinates.tolist(),
                    SCF_converged=True, surface_points=len(mf.with_solvent.surface["charge_exp"]))
                atomic_write_json(raw / "baseline.json", baseline)
                response = pcm.analytical_grad_vmat(mf.with_solvent, dm, mo, occupied)
                np.save(raw / "analytic_response_MO_occ.npy", response.get())
                reference = {(a,x): response[a,x].get() for a,x in probes}
                to_cpu = lambda value: value.get() if hasattr(value, "get") else np.asarray(value)
                host_dm, host_mo, host_occ = to_cpu(dm), to_cpu(mo), to_cpu(occupied)
                del response, mf, dm, mo, occupied
                gc.collect()
                dft._free_gpu_pool()
                checks = []
                for atom, axis in probes:
                    for step in plan["settings"]["steps_Bohr"]:
                        values = []
                        for sign in (-1,1):
                            shifted = coordinates.copy()
                            shifted[atom,axis] += sign*step
                            mf = mean_field(shifted)
                            _, matrix = mf.with_solvent._get_vind(cupy.asarray(host_dm))
                            projection = cupy.asarray(host_mo).T @ matrix @ cupy.asarray(host_mo[:, host_occ > 0])
                            value = projection.get()
                            name = f"atom{atom}_axis{axis}_h{step}_sign{sign}"
                            np.save(raw / (name+".npy"), value)
                            atomic_write_json(raw / (name+".json"), dict(coordinates_Bohr=shifted.tolist(),
                                fixed_AO_density=True, fixed_MO_coefficients=True, SCF_performed=False))
                            values.append(value)
                            del projection, matrix, mf
                            gc.collect()
                            dft._free_gpu_pool()
                        numeric = (values[1]-values[0])/(2*step)
                        analytic = reference[atom,axis]
                        error = float(np.linalg.norm(numeric-analytic)/max(1e-14,np.linalg.norm(analytic)))
                        checks.append(dict(atom=atom, axis=axis, step_Bohr=step, relative_error=error,
                            analytic_norm=float(np.linalg.norm(analytic)), numeric_norm=float(np.linalg.norm(numeric))))
                        atomic_write_json(raw / "response_checks.json", checks)
                        status["completed_response_checks"] = len(checks)
                        atomic_write_json(status_path, status)
                if any(not np.isfinite(r["relative_error"]) or r["relative_error"] > .01 for r in checks):
                    raise ValueError("Independent fixed-density response check failed")
        if monitor.peak > plan["acceptance"]["maximum_GPU_MiB"]:
            raise ValueError("Measured resource footprint exceeds the declared device profile")
        api.verify_sources(plan["producer_sha256"])
        api.verify_sources(plan["input_sha256"])
        status.update(status="response_resource_and_sampled_derivatives_passed", checks=checks,
            output_sha256={str(p): api.digest(p) for p in raw.rglob("*") if p.is_file()})
    except BaseException as error:
        status.update(status="response_pilot_failed_preserved", error=repr(error))
        raise
    finally:
        status.update(finished=time.time(), observed_peak_GPU_MiB=monitor.peak if monitor else None)
        atomic_write_json(status_path, status)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("action", choices=("prepare", "run"))
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--source", type=Path)
    p.add_argument("--controls", type=Path)
    args = p.parse_args()
    if args.action == "prepare":
        if not args.source or not args.controls:
            p.error("Preparation requires the checked parent and full spin controls")
        from redox.workflows.resources import allocation_guard, host_capacity, cpu_allocation
        with allocation_guard():
            host_capacity(1)
            os.sched_setaffinity(0,cpu_allocation(1))
        prepare(args)
    else:
        run(args)


if __name__ == "__main__":
    main()
