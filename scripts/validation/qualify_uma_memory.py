#!/usr/bin/env python
"""Measure a bounded UMA preoptimization workload inside a system GPU reservation.

This qualification concerns UMA only. It does not admit DFT or Hessians to small
devices. The same script can execute another in-scope motif after qualification;
its per-process allocation cap and input fingerprints remain part of the evidence.
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import threading
import time

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from redox.core.common import read_manifest, LIBRARY, UMA
from redox.core.io import atomic_write_json
from redox.workflows.resources import verify_gate_owner


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--only", required=True, help="registered motif or motif:state")
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    gid, _, selected_state = args.only.partition(":")
    verify_gate_owner()
    destination = args.out.resolve(); destination.mkdir(parents=True, exist_ok=False)
    report = dict(status="initializing", workload="UMA gas preoptimization only", id=gid, selected_state=selected_state or None,
                  scientific_model_changed=False, final_solution_properties_validated=False,
                  pid=os.getpid(), process_start=Path("/proc/self/stat").read_text().rsplit(")", 1)[1].split()[19],
                  cpu_affinity=sorted(os.sched_getaffinity(0)), started=time.time(),
                  allocation_fraction=.75, states=[])
    path = destination / "status.json"
    atomic_write_json(path, report)
    stop = threading.Event(); sampler = None; samples = []
    try:
        import torch
        from redox.qm import uma
        rows = [row for row in read_manifest() if row["id"] == gid and (not selected_state or row["state"] == selected_state)]
        if not rows:
            raise ValueError("No registered states selected")
        gpu = int(os.environ["CUDA_VISIBLE_DEVICES"])
        if torch.cuda.device_count() != 1:
            raise RuntimeError("Expected one system-reserved GPU")
        device = torch.cuda.get_device_properties(0)
        if not 14*1024**3 <= device.total_memory <= 17*1024**3:
            raise RuntimeError("This study qualifies the 16-GB class explicitly")
        torch.cuda.set_per_process_memory_fraction(.75, 0)
        torch.cuda.reset_peak_memory_stats(0)
        confs = sorted((LIBRARY / gid / "conformers").glob("conf_*.xyz"))
        if not confs:
            raise RuntimeError("No prepared conformer ensemble")
        identity = uma.model_identity(uma.DEFAULT_MODEL)
        report.update(status="computing", torch_version=torch.__version__, model_identity=identity,
                      producer_sha256=uma.file_sha256(Path(uma.__file__)),
                      runtime_sha256=uma.file_sha256(Path(__file__)), physical_gpu=gpu,
                      device_name=device.name, total_memory_bytes=device.total_memory,
                      compute_capability=[device.major, device.minor],
                      max_atoms=max(int(conf.read_text().splitlines()[0]) for conf in confs),
                      seed_identity=uma.seed_identity(confs))
        atomic_write_json(path, report)

        def sample():
            while not stop.is_set():
                try:
                    text = subprocess.check_output(
                        ["nvidia-smi", f"--id={gpu}", "--query-gpu=memory.used",
                         "--format=csv,noheader,nounits"], text=True, timeout=5)
                    samples.append(dict(time=time.time(), used_MiB=float(text.strip())))
                except (OSError, ValueError, subprocess.SubprocessError):
                    pass
                stop.wait(.5)

        sampler = threading.Thread(target=sample, daemon=True); sampler.start()
        previous = sys.argv
        try:
            sys.argv = ["redox.qm.uma", "--only", args.only, "--device", "cuda"]
            uma.main()
        finally:
            sys.argv = previous
        for row in rows:
            directory = UMA / row["id"] / row["state"]
            result = json.loads((directory / "result.json").read_text())
            errors = uma.cache_errors(result, row, directory / "relaxed.xyz", confs, identity, .02)
            if errors:
                raise RuntimeError(f"Incomplete/unverifiable UMA result: {row['state']} {errors}")
            report["states"].append(dict(state=row["state"], chosen_mult=result["chosen_mult"],
                                         energy_eV=result["energy_eV"], converged=result["converged"],
                                         result_sha256=uma.file_sha256(directory / "result.json")))
        torch.cuda.synchronize()
        report.update(peak_allocated_bytes=torch.cuda.max_memory_allocated(0),
                      peak_reserved_bytes=torch.cuda.max_memory_reserved(0),
                      peak_device_used_MiB=max((entry["used_MiB"] for entry in samples), default=None))
        if not samples or report["peak_device_used_MiB"]*1024**2 > .8*device.total_memory:
            raise RuntimeError("UMA qualification lacks measured memory headroom")
        report["status"] = "bounded_UMA_preoptimization_verified"
    except BaseException as exc:
        report.update(status="failed_preserving_outputs", error=repr(exc))
        raise
    finally:
        stop.set()
        if sampler is not None:
            sampler.join(timeout=6)
        report.setdefault("peak_device_used_MiB", max((entry["used_MiB"] for entry in samples), default=None))
        report.update(finished=time.time(), memory_samples=samples)
        atomic_write_json(path, report)


if __name__ == "__main__":
    main()
