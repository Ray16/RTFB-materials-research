#!/usr/bin/env python
"""Finish the last frozen parent benchmark's source-addressed spin check."""

import argparse
from contextlib import contextmanager
import json
import os
from pathlib import Path
import socket
import time

from redox.core.common import ROOT, DFT
from redox.core.io import atomic_write_json
from redox.workflows.chunked import file_hash

GID = "bm_n_n_bis_2_2_methoxyethoxy_ethyl_pyromellitic_diimide"
STATE = "red2"
OUT = ROOT / "results/validation/publication_benchmark/spin_completion_2026-10-21"
SCRIPT = Path(__file__).resolve()


def prepare():
    from redox.validation.spins import snapshot, declared
    from redox.validation.spin_evidence import state_spin_evidence, spin_comparison_running
    from redox.qm.native_thermal import load_native_thermal
    from redox.core.protocol import record_path

    if spin_comparison_running(GID):
        raise RuntimeError("An existing spin writer is still active")
    source = snapshot(GID, STATE)
    if source is None or state_spin_evidence(GID, STATE)["spin_confidence"] != "DFT_check_required":
        raise ValueError("This continuation requires the remaining exact benchmark spin hold")
    state = DFT / GID / STATE
    if load_native_thermal(state, source["charge"], source["mult"]) is None:
        raise ValueError("Primary native minimum must qualify first")
    OUT.mkdir()
    frozen = OUT / "chosen_snapshot"
    frozen.mkdir()
    for label, path in [("opt.xyz", state / "opt.xyz"), ("result.json", state / "result.json"),
                        ("active_SP.json", record_path(state, source["charge"], source["mult"]))]:
        (frozen / label).write_bytes(path.read_bytes())
    atomic_write_json(frozen / "manifest_row.json", declared(GID, STATE))
    paths = [SCRIPT] + [ROOT / "src/redox" / name for name in [
        "validation/spins.py", "validation/spin_evidence.py", "qm/dft.py", "qm/sp.py",
        "qm/native_thermal.py", "qm/projected_thermal.py", "qm/qrrho.py", "qm/native_thermal_repair.py",
        "workflows/geometry_checkpoints.py", "core/common.py", "core/protocol.py"]]
    plan = dict(role="Retrospective exact-parent method control; no synthesis admission", id=GID, state=STATE,
                chosen_identity=source, producer_sha256={str(p): file_hash(p) for p in paths},
                input_sha256={str(p): file_hash(p) for p in frozen.iterdir()},
                numerical_settings="Unchanged source-addressed alternate-spin implementation",
                durable_geometry_checkpoints=True, native_gates_unchanged=True)
    atomic_write_json(OUT / "plan.json", plan)


@contextmanager
def checkpoint_optimizers():
    from redox.qm import dft
    from redox.workflows.geometry_checkpoints import checkpoint_callback, write_geometry

    previous = dft._geomopt_fn
    index = 0

    def driver(backend):
        kernel = previous(backend)

        def run(method, **kwargs):
            nonlocal index
            folder = OUT / "geometry_checkpoints" / f"optimizer_{index}"
            index += 1
            folder.mkdir(parents=True)
            save = checkpoint_callback(folder)
            original = kwargs.get("callback")

            def callback(environment):
                save(environment)
                if original:
                    original(environment)

            kwargs["callback"] = callback
            result = kernel(method, **kwargs)
            write_geometry(result[1], folder / "optimizer_returned.xyz", "Unqualified optimizer return")
            atomic_write_json(folder / "optimizer_returned.json", dict(optimizer_converged=bool(result[0]),
                              physical_admission="Requires normal spin, topology, SP and native gates"))
            return result

        return run

    dft._geomopt_fn = driver
    try:
        yield
    finally:
        dft._geomopt_fn = previous


def worker():
    from redox.workflows.resources import verify_gate_owner
    from redox.validation.spins import snapshot, check
    from redox.validation.spin_evidence import state_spin_evidence

    verify_gate_owner()
    if len(os.sched_getaffinity(0)) != 2:
        raise RuntimeError("Two disjoint CPU slots are required")
    caps = ("OMP_NUM_THREADS", "OMP_THREAD_LIMIT", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS", "XTB_THREADS")
    if any(os.environ.get(key) != "2" for key in caps):
        raise RuntimeError("All native pools must be capped at two")
    plan = json.loads((OUT / "plan.json").read_text())
    for group in ("producer_sha256", "input_sha256"):
        if any(file_hash(p) != expected for p, expected in plan[group].items()):
            raise ValueError("Spin continuation source changed")
    if snapshot(GID, STATE) != plan["chosen_identity"]:
        raise ValueError("Chosen physical state changed before the continuation")
    status = dict(status="computing", owner=dict(host=socket.gethostname(), pid=os.getpid(),
                  process_start=Path("/proc/self/stat").read_text().rsplit(")", 1)[1].split()[19],
                  gpu=os.environ["CUDA_VISIBLE_DEVICES"], cpus=sorted(os.sched_getaffinity(0))), started=time.time())
    atomic_write_json(OUT / "status.json", status)
    try:
        with checkpoint_optimizers():
            result = check(GID, STATE)
        evidence = state_spin_evidence(GID, STATE)
        status.update(result=result, publication_spin_evidence=evidence,
                      status="spin_admission_verified" if evidence["spin_primary_admissible"] else "spin_check_requires_review")
        if result["status"] in ("failed", "source_stale"):
            raise RuntimeError("Alternate-spin calculation failed; outputs preserved")
    except BaseException as error:
        status.update(status="failed_preserved", error=repr(error))
        raise
    finally:
        status["finished"] = time.time()
        atomic_write_json(OUT / "status.json", status)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("prepare", "worker"))
    args = parser.parse_args()
    prepare() if args.action == "prepare" else worker()
