#!/usr/bin/env python
"""Validate, graft, compute and audit registered motifs. Default is a read-only plan.

  python scripts/pipeline/run_motifs.py --config config/custom_candidates.py --out results/run_name
  python scripts/pipeline/run_motifs.py --config config/custom_candidates.py --out results/run_name --execute --gpu 0

GPU steps always run through gpu_reserve. Existing outputs are reused only when valid;
repair that would replace an existing geometry must be reviewed separately.
"""
import argparse
import fcntl
import hashlib
import json
import os
from pathlib import Path
import shlex
import shutil
import socket
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from redox.build.candidates import _load_config
from redox.build.contracts import validate_group
from redox.core.io import atomic_write_json
from redox.workflows.resources import (cpu_allocation, GATE, allocation_guard,
                                       gpu_suitability, host_capacity, free_slots)
from redox.screening.scorecard import CANDIDATE_CONFIGS


def selection_hash(scaffold, groups):
    payload = dict(scaffold=scaffold, groups=sorted(groups, key=lambda group: group["id"]))
    return hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest()


def resume_config_matches(previous, config, scaffold, groups):
    """Permit append-only registry growth while binding the selected definitions.

    Older runs can be checked against an exact archived registry made by the
    append-only registration workflow. A changed selected molecule still fails.
    """
    selected = selection_hash(scaffold, groups)
    if previous.get("selected_definitions_sha256"):
        return previous["selected_definitions_sha256"] == selected
    old_hash = previous.get("input_sha256")
    if old_hash == hashlib.sha256(config.read_bytes()).hexdigest():
        return True
    if not old_hash or len(old_hash) != 64 or any(c not in "0123456789abcdef" for c in old_hash):
        return False
    ids = {group["id"] for group in groups}
    for archive in (ROOT / "results/discovery").glob(f"*/config_before/{old_hash}.py"):
        if hashlib.sha256(archive.read_bytes()).hexdigest() != old_hash:
            continue
        old_scaffold, old_groups = _load_config(archive)
        selected_old = [group for group in old_groups if group["id"] in ids]
        if {group["id"] for group in selected_old} == ids and selection_hash(old_scaffold, selected_old) == selected:
            return True
    return False


def run_stage(argv, log, env, retry_delays=(1, 2, 4, 8, 16)):
    """Retry only a gate refusal that proves the scientific command never ran.

    A recently exited CUDA process can leave the driver reporting memory briefly.
    Every retry still goes through the atomic gate; real stage failures stop at once.
    """
    attempts = []
    for attempt in range(len(retry_delays) + 1):
        attempt_log = log if attempt == 0 else log.with_name(
            f"{log.stem}_gate_retry_{attempt}{log.suffix}")
        started = time.time()
        with attempt_log.open("w") as stream:
            proc = subprocess.run(argv, cwd=ROOT, env=env, stdout=stream, stderr=subprocess.STDOUT)
        refused = (argv[:2] == [str(GATE), "run"] and proc.returncode == 1
                   and attempt_log.read_text().strip()
                   == f"gpu_reserve: gpu{argv[2]} BUSY/reserved -> not running")
        attempts.append(dict(log=str(attempt_log), started=started, finished=time.time(),
                             returncode=proc.returncode, reservation_refused=refused))
        if not refused or attempt == len(retry_delays):
            return proc, attempts
        print(f"GPU reservation unavailable; retrying in {retry_delays[attempt]}s", flush=True)
        time.sleep(retry_delays[attempt])



def commands(config, ids, output, backend, threads, nconf, thermal_definition="native"):
    py = sys.executable
    jobs = []
    for gid in ids:
        frequency_args = ["--no-freq"] if thermal_definition == "native" else []
        jobs.extend([
            ("build", False, [py, "-m", "redox.build.candidates", "--config", str(config), "--only", gid]),
            ("uma", True, [py, "-m", "redox.qm.uma", "--only", gid, "--device", "cuda" if backend == "gpu" else "cpu"]),
            ("dft", True, [py, "-m", "redox.qm.dft", "--only", gid, "--backend", backend,
                           "--nthreads", str(threads), "--nconf", str(nconf), *frequency_args]),
            ("sp", True, [py, "-m", "redox.qm.sp", "--only", gid, "--backend", backend]),
        ])
        if thermal_definition == "native":
            jobs.append(("native_thermal", True, [py, "-m", "redox.qm.native_thermal", "--only",
                                                  f"ferrocene,{gid}", "--threads", str(threads),
                                                  "--select-policy", "--repair-minimum"]))
            jobs.append(("fc_reference_sync", False, [py, str(ROOT / "scripts/pipeline/set_fc_reference.py")]))
        jobs.extend([
            ("coverage", False, [py, "-m", "redox.validation.readiness", "--only", gid, "--strict"]),
            ("cross_points", True, [py, "-m", "redox.properties.reorg", "--only", gid, "--backend", backend]),
        ])
    jobs.append(("review_tables", False, [py, "-m", "redox.workflows.tables",
                                          "--out", str(output / "tables")]))
    return jobs


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", type=Path, default=ROOT / "config/custom_candidates.py")
    ap.add_argument("--only", help="comma-separated registered IDs")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--execute", action="store_true")
    ap.add_argument("--resume", action="store_true", help="resume a preserved run with matching inputs")
    ap.add_argument("--skip-review", action="store_true", help="fleet mode: verify selected records; collect tables separately")
    ap.add_argument("--backend", choices=["gpu", "cpu"], default="gpu")
    ap.add_argument("--gpu", type=int)
    ap.add_argument("--threads", type=int, default=2)
    ap.add_argument("--nconf", type=int, default=3)
    ap.add_argument("--cpus", help="comma-separated CPU indices; otherwise choose unallocated CPUs")
    ap.add_argument("--thermal-definition", choices=["native", "legacy"], default="native",
                    help="new computations use derivative-checked solution DFT thermochemistry")
    a = ap.parse_args()
    if a.resume and not a.execute:
        ap.error("--resume requires --execute")
    config, output = a.config.resolve(), a.out.resolve()
    allowed = {ROOT / "config" / f"{name}.py" for _, name in CANDIDATE_CONFIGS}
    allowed.add(ROOT / "config/redox_groups.py")
    if config not in allowed:
        ap.error("use a registered candidate config (custom_candidates.py for new motifs)")
    scaffold, groups = _load_config(config)
    if a.only:
        ids = set(a.only.split(","))
        unknown = ids - {g["id"] for g in groups}
        if unknown:
            ap.error(f"unknown candidate IDs: {sorted(unknown)}")
        groups = [g for g in groups if g["id"] in ids]
    if not groups:
        ap.error("config contains no selected motifs; add GROUPS entries first")
    if a.threads < 1 or a.nconf < 1:
        ap.error("threads and nconf must be positive")
    for g in groups:
        validate_group(scaffold, g)
    ids = [g["id"] for g in groups]
    if a.backend == "cpu" and a.thermal_definition == "native":
        ap.error("native thermochemistry is validated for GPU execution; select a GPU backend")
    jobs = commands(config, ids, output, a.backend, a.threads, a.nconf, a.thermal_definition)
    if a.skip_review:
        jobs = [job for job in jobs if job[0] not in ("review_tables", "fc_reference_sync")]
    if not a.execute:
        for label, gpu, argv in jobs:
            prefix = f"{GATE} run <gpu> -- " if gpu and a.backend == "gpu" else ""
            print(f"{label}: {prefix}{shlex.join(argv)}")
        print("Plan only. Add --execute and --gpu INDEX to run with a reservation.")
        return
    previous = None
    if output.exists():
        if not a.resume:
            ap.error("choose a new output directory or use --resume; existing runs are preserved")
        try:
            previous = json.loads((output / "run_manifest.json").read_text())
        except (OSError, ValueError):
            ap.error("resume requires an intact prior run manifest")
        expected = dict(config=str(config), ids=ids, backend=a.backend, nconf=a.nconf,
                        thermal_definition=a.thermal_definition)
        if (any(previous.get(key) != value for key, value in expected.items())
                or not resume_config_matches(previous, config, scaffold, groups)):
            ap.error("resume inputs or methods differ from the preserved run")
    elif a.resume:
        ap.error("there is no preserved run to resume")
    if a.backend == "gpu" and (a.gpu is None or not GATE.is_file()):
        ap.error("GPU execution requires --gpu and the system reservation gate")
    allowed_cpus = sorted(os.sched_getaffinity(0))
    if (a.threads > min(len(allowed_cpus), max(1, (os.cpu_count() or 1) // 2))
            or os.getloadavg()[0] + a.threads > (os.cpu_count() or 1)):
        ap.error("insufficient CPU capacity on this shared node")
    taskset = shutil.which("taskset")
    if not taskset:
        ap.error("taskset is required to constrain worker CPU affinity")
    try:
        requested = [int(v) for v in a.cpus.split(",")] if a.cpus else None
        with allocation_guard():
            cpu_set = cpu_allocation(a.threads, requested)
            node_preflight = host_capacity(a.threads)
            if a.backend == "gpu":
                node_preflight["gpu"] = gpu_suitability(a.gpu)
                host = socket.gethostname().split(".")[0]
                if (host, a.gpu) not in free_slots([host]):
                    raise RuntimeError("GPU occupied or reserved; no candidate calculation started")
            original_affinity = set(os.sched_getaffinity(0))
            os.sched_setaffinity(0, cpu_set)
    except (ValueError, RuntimeError, OSError, subprocess.SubprocessError) as exc:
        ap.error(str(exc))
    env = dict(os.environ, PYTHONPATH=str(ROOT / "src"))
    env["REDOX_PIPELINE_OWNER_PID"] = str(os.getpid())
    env["REDOX_PIPELINE_OWNER_START"] = Path("/proc/self/stat").read_text().rsplit(")", 1)[1].split()[19]
    env["PATH"] = str(Path(sys.executable).parent) + os.pathsep + env.get("PATH", "")
    xtb = env.get("XTB_BIN") or shutil.which("xtb", path=env["PATH"])
    if not xtb:
        ap.error("the selected environment has no xTB executable")
    env["XTB_BIN"] = xtb
    for key in ("OMP_NUM_THREADS", "OMP_THREAD_LIMIT", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
        env[key] = str(a.threads)
    import site
    libs = [str(p) for base in site.getsitepackages() for p in Path(base).glob("nvidia/*/lib")]
    env["LD_LIBRARY_PATH"] = ":".join(libs + [env.get("LD_LIBRARY_PATH", "")])
    # One writer owns each motif for the entire run; kernel releases locks on crash/exit.
    handles = []
    try:
        lockdir = ROOT / "logs/pipeline_locks"; lockdir.mkdir(parents=True, exist_ok=True)
        for gid in sorted(ids):
            f = (lockdir / f"{gid}.lock").open("a+")
            handles.append(f)
            try:
                fcntl.flock(f, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError:
                ap.error(f"{gid}: another pipeline run owns this motif")
            f.seek(0); f.truncate()
            f.write(json.dumps(dict(host=socket.gethostname(), pid=os.getpid(), id=gid)))
            f.flush()
        output.mkdir(parents=True, exist_ok=a.resume)
        attempt_output = output
        if previous is not None:
            attempt_output = output / f"resume_{time.time_ns()}"
            attempt_output.mkdir()
            shutil.copy2(output / "run_manifest.json", attempt_output / "previous_run_manifest.json")
            jobs = commands(config, ids, attempt_output, a.backend, a.threads, a.nconf, a.thermal_definition)
            if a.skip_review:
                jobs = [job for job in jobs if job[0] not in ("review_tables", "fc_reference_sync")]
        manifest = dict(config=str(config), input_sha256=hashlib.sha256(config.read_bytes()).hexdigest(),
                        selected_definitions_sha256=selection_hash(scaffold, groups),
                        ids=ids, host=socket.gethostname(), pid=os.getpid(), python=sys.executable,
                        backend=a.backend, gpu=a.gpu, threads=a.threads, cpu_affinity=cpu_set,
                        node_preflight=node_preflight,
                        xtb=xtb, nconf=a.nconf, thermal_definition=a.thermal_definition,
                        prior_run_manifest=str(attempt_output / "previous_run_manifest.json") if previous else None,
                        attempt_output=str(attempt_output), steps=[])
        mp = output / "run_manifest.json"
        for index, (label, compute, argv) in enumerate(jobs):
            argv = [taskset, "--cpu-list", ",".join(map(str, cpu_set)), *argv]
            if compute and a.backend == "gpu":
                argv = [str(GATE), "run", str(a.gpu), "--", *argv]
            log = attempt_output / f"{index:02d}_{label}.log"
            entry = dict(stage=label, argv=argv, log=str(log), started=time.time())
            manifest["steps"].append(entry)
            atomic_write_json(mp, manifest)
            print(f"[{index + 1}/{len(jobs)}] {label}: {log}", flush=True)
            proc, attempts = run_stage(argv, log, env)
            entry.update(returncode=proc.returncode, finished=time.time(), attempts=attempts)
            atomic_write_json(mp, manifest)
            if proc.returncode:
                from redox.validation.readiness import inspect_states
                from redox.core.common import read_manifest
                manifest.update(status="resource_unavailable" if attempts[-1]["reservation_refused"] else "failed",
                                failed_stage=label,
                                scientific_accuracy_validated=False,
                                selected_state_qc=inspect_states([r for r in read_manifest() if r["id"] in ids]))
                atomic_write_json(mp, manifest)
                raise SystemExit(f"{label} failed; see {log}")
        # Written last, after the independent audit passes; QC limits remain in tables.
        from redox.validation.readiness import inspect_states
        from redox.core.common import read_manifest
        states = inspect_states([r for r in read_manifest() if r["id"] in ids])
        if a.skip_review:
            from redox.properties.potentials import fc_live, thermal_definition
            from redox.core.common import read_result, free_energy
            reference=fc_live()
            definition=thermal_definition(read_result("ferrocene","ox"))
            selected=[r for r in read_manifest() if r["id"] in ids]
            records=[read_result(r["id"],r["state"]) for r in selected]
            if (reference is None or any(s["status"]!="ready" for s in states)
                    or any(free_energy(r)[1]!="ok" or free_energy(r)[0] is None
                           or thermal_definition(r)!=definition for r in records)):
                raise SystemExit("Selected records or matched Fc reference remain unvalidated")
            manifest.update(status="selected_outputs_verified",selected_state_qc=states,
                            selected_states_ready=True,scientific_accuracy_validated=False,
                            table_collection_pending=True)
            atomic_write_json(mp, manifest)
            print(f"Verified selected records: {mp}; tables await separate collection")
            return
        review = json.loads((attempt_output / "tables/review_status.json").read_text())
        if review.get("arithmetic_discrepancies") != 0 or not (attempt_output / "tables/Candidates_computed.xlsx").is_file():
            raise SystemExit("Final property artifacts are missing or failed the independent audit")
        import csv
        with (attempt_output / "tables/redox_potentials.csv").open() as stream:
            selected_events = [r for r in csv.DictReader(stream) if r["id"] in ids]
        if ({r["id"] for r in selected_events} != set(ids)
                or any(r.get("status") != "ok" for r in selected_events)
                or any(r["status"] != "ready" for r in states)):
            raise SystemExit("Selected motifs still have incomplete or unvalidated property records")
        manifest.update(status="tables_audited", selected_state_qc=states,
                        selected_states_ready=all(r["status"] == "ready" for r in states),
                        scientific_accuracy_validated=False)
        atomic_write_json(mp, manifest)
        print(f"Audited property tables: {attempt_output / 'tables'}")
    finally:
        for f in handles:
            f.close()
        os.sched_setaffinity(0, original_affinity)


if __name__ == "__main__":
    main()
