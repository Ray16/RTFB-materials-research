#!/usr/bin/env python
"""Validate, graft, compute and audit registered motifs. Default is a read-only plan.

  python scripts/pipeline/run_motifs.py --config config/custom_candidates.py --out results/run_name
  python scripts/pipeline/run_motifs.py --config config/custom_candidates.py --out results/run_name --execute --gpu 0
  ... --execute --gpu 0,1,2   # optimize the redox states of each motif concurrently, one GPU each

GPU steps always run through gpu_reserve. Existing outputs are reused only when valid;
repair that would replace an existing geometry must be reviewed separately.
"""
import argparse
from dataclasses import dataclass
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

from redox.core.common import ROOT
from redox.structures.candidates import _load_config
from redox.structures.contracts import validate_group
from redox.core.io import atomic_write_json
from redox.execution.resources import (cpu_allocation, GATE, allocation_guard,
                                       gpu_suitability, host_capacity, free_slots)
from redox.screening.scorecard import CANDIDATE_CONFIGS
from redox.workflows.stages import execute_graph, validate_graph


def selection_hash(scaffold, groups):
    payload = dict(scaffold=scaffold, groups=sorted(groups, key=lambda group: group["id"]))
    return hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest()


def resume_config_matches(previous, config, scaffold, groups, *, root=None):
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
    for archive in ((ROOT if root is None else root) / "results/discovery").glob(f"*/config_before/{old_hash}.py"):
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


def state_jobs(argv, gid, states):
    """One command per redox state: the DFT optimizations of a motif's states are independent
    (each writes only calcs/dft/<id>/<state>/), so they can run on separate GPUs."""
    at = argv.index("--only") + 1
    if argv[at] != gid:
        raise ValueError(f"expected --only {gid}")
    return [(state, [*argv[:at], f"{gid}:{state}", *argv[at + 1:]]) for state in states]


@dataclass(frozen=True)
class MotifNode:
    """One unit of work in the campaign DAG scheduled by stages.execute_graph.

    name/dependencies are what execute_graph consumes; label/argv/compute reuse the
    exact per-stage command from commands() verbatim (byte-identical argv -> zero
    science change). gid/state/states and lock carry the metadata the resume
    predicate and the one-time ferrocene prerequisite need.
    """
    name: str
    label: str
    argv: tuple
    compute: bool
    dependencies: tuple = ()
    gid: str = ""
    state: str = ""
    states: tuple = ()
    lock: str = ""


def _split_by_label(per_candidate_jobs):
    by_label = {}
    for label, compute, argv in per_candidate_jobs:
        by_label.setdefault(label, []).append((compute, list(argv)))
    return by_label


def build_graph(config, groups, output, backend, threads, nconf,
                thermal_definition="native", ferrocene_states=()):
    """Bag-of-tasks DAG for the whole pool, one node per unit of GPU/CPU work.

    Replaces the flat per-candidate serial job list. Per-(candidate, redox-state)
    DFT nodes let independent states and candidates pipeline across GPU slots
    instead of the slowest state (ox1) gating the next candidate. Edges encode the
    scientific order exactly: uma -> dft(per state) -> sp -> native_thermal and
    sp -> coverage -> cross_points (coverage gates reorg). Ferrocene is hoisted to
    a single prerequisite node (its own lock) and fc_reference_sync is one node after
    every native_thermal, so the shared calcs/dft/ferrocene and config/project.json
    are never written concurrently under cross-candidate parallelism.

    Review-table building is deliberately NOT a node here; it is taken off the
    campaign critical path and run as a separate downstream job.
    """
    native = thermal_definition == "native"
    nodes = []
    native_thermal_names = []
    fc_argv = None
    fc_compute = False
    for g in groups:
        gid = g["id"]
        states = tuple(s[0] for s in g["states"])
        by_label = _split_by_label(
            commands(config, [gid], output, backend, threads, nconf, thermal_definition))
        b_compute, b_argv = by_label["build"][0]
        nodes.append(MotifNode(f"build:{gid}", "build", tuple(b_argv), b_compute,
                               (), gid=gid))
        u_compute, u_argv = by_label["uma"][0]
        nodes.append(MotifNode(f"uma:{gid}", "uma", tuple(u_argv), u_compute,
                               (f"build:{gid}",), gid=gid, states=states))
        d_compute, d_argv = by_label["dft"][0]
        dft_names = []
        for state, sargv in state_jobs(d_argv, gid, states):
            name = f"dft:{gid}:{state}"
            dft_names.append(name)
            nodes.append(MotifNode(name, "dft", tuple(sargv), d_compute,
                                   (f"uma:{gid}",), gid=gid, state=state))
        s_compute, s_argv = by_label["sp"][0]
        nodes.append(MotifNode(f"sp:{gid}", "sp", tuple(s_argv), s_compute,
                               tuple(dft_names), gid=gid, states=states))
        coverage_deps = [f"sp:{gid}"]
        if native:
            nt_compute, nt_argv = by_label["native_thermal"][0]
            at = nt_argv.index("--only") + 1
            # commands() targets "ferrocene,<gid>"; split so the shared ferrocene
            # namespace is computed exactly once by its own prerequisite node and
            # each candidate's native_thermal only touches its own records.
            if fc_argv is None:
                fc_compute, fc_argv = by_label["fc_reference_sync"][0]
                ferrocene_argv = [*nt_argv[:at], "ferrocene", *nt_argv[at + 1:]]
                nodes.append(MotifNode("native_thermal:ferrocene", "native_thermal",
                                       tuple(ferrocene_argv), nt_compute, (),
                                       gid="ferrocene", states=tuple(ferrocene_states),
                                       lock="ferrocene"))
            cand_argv = [*nt_argv[:at], gid, *nt_argv[at + 1:]]
            nt_name = f"native_thermal:{gid}"
            nodes.append(MotifNode(nt_name, "native_thermal", tuple(cand_argv), nt_compute,
                                   (f"sp:{gid}", "native_thermal:ferrocene"),
                                   gid=gid, states=states))
            native_thermal_names.append(nt_name)
            # readiness --strict requires the candidate's native thermal record.
            coverage_deps.append(nt_name)
        c_compute, c_argv = by_label["coverage"][0]
        nodes.append(MotifNode(f"coverage:{gid}", "coverage", tuple(c_argv), c_compute,
                               tuple(coverage_deps), gid=gid))
        r_compute, r_argv = by_label["cross_points"][0]
        nodes.append(MotifNode(f"cross_points:{gid}", "cross_points", tuple(r_argv), r_compute,
                               (f"coverage:{gid}", f"sp:{gid}"), gid=gid))
    if native and fc_argv is not None:
        # One reference sync after EVERY native_thermal (ferrocene included), before
        # anything downstream reads the Fc reference (final verify / separate tables).
        nodes.append(MotifNode("fc_sync", "fc_reference_sync", tuple(fc_argv), fc_compute,
                               tuple(["native_thermal:ferrocene", *native_thermal_names])))
    validate_graph(nodes)
    return nodes


def node_done(node):
    """True only when a node's scientific output already exists AND is valid, so a
    resume skips it with no subprocess and no GPU reservation (resume is O(remaining)).

    Conservative by construction: anything uncertain returns False and is simply
    re-run; the stage subprocesses themselves reuse valid outputs idempotently, so a
    False never corrupts science, while a True must be certain. Only the expensive
    GPU stages with a cheap unambiguous completion check are reported done.
    """
    from redox.core.common import DFT, UMA, read_result, free_energy
    gid, states = node.gid, node.states
    if node.label == "dft":
        sd = DFT / gid / node.state
        if not ((sd / "result.json").is_file() and (sd / "opt.xyz").is_file()):
            return False
        try:
            from redox.backends.dft import _opt_protocol_mismatch
            r = json.loads((sd / "result.json").read_text())
            return not _opt_protocol_mismatch(r, int(r["charge"]))
        except (OSError, ValueError, KeyError):
            return False
    if node.label == "uma":
        return bool(states) and all((UMA / gid / s / "result.json").is_file() for s in states)
    if node.label == "sp":
        if not states:
            return False
        recs = [read_result(gid, s) for s in states]
        return all(r is not None and r.get("sp_status") == "ok" for r in recs)
    if node.label == "native_thermal":
        if not states:
            return False
        for s in states:
            r = read_result(gid, s)
            if (r is None or r.get("thermal_definition") != "native_solution_dft_qrrho_spin"
                    or r.get("thermal_qc") != "ok"):
                return False
        return True
    return False


def commands(config, ids, output, backend, threads, nconf, thermal_definition="native"):
    py = sys.executable
    jobs = []
    for gid in ids:
        frequency_args = ["--no-freq"] if thermal_definition == "native" else []
        jobs.extend([
            ("build", False, [py, "-m", "redox.structures.candidates", "--config", str(config), "--only", gid]),
            ("uma", True, [py, "-m", "redox.backends.uma_current", "--stage", "uma", "--only", gid, "--device", "cuda" if backend == "gpu" else "cpu"]),
            ("dft", True, [py, "-m", "redox.backends.uma_current", "--stage", "dft", "--only", gid, "--backend", backend,
                           "--nthreads", str(threads), "--nconf", str(nconf), *frequency_args]),
            ("sp", True, [py, "-m", "redox.backends.sp", "--only", gid, "--backend", backend]),
        ])
        if thermal_definition == "native":
            jobs.append(("native_thermal", True, [py, "-m", "redox.thermochemistry.native", "--only",
                                                  f"ferrocene,{gid}", "--threads", str(threads),
                                                  "--select-policy", "--repair-minimum"]))
            jobs.append(("fc_reference_sync", False, [py, str(ROOT / "scripts/pipeline/set_fc_reference.py")]))
        jobs.extend([
            ("coverage", False, [py, "-m", "redox.quality.readiness", "--only", gid, "--strict"]),
            ("cross_points", True, [py, "-m", "redox.properties.reorg", "--only", gid, "--backend", backend]),
        ])
    jobs.append(("review_tables", False, [py, "-m", "redox.reporting.tables",
                                          "--out", str(output / "tables")]))
    return jobs


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", type=Path, default=ROOT / "config/custom_candidates.py")
    ap.add_argument("--only", help="comma-separated registered IDs")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--execute", action="store_true")
    ap.add_argument("--resume", action="store_true", help="resume a preserved run with matching inputs")
    ap.add_argument("--skip-review", action="store_true", help="fleet mode: verify selected records; collect tables separately")
    ap.add_argument("--backend", choices=["gpu", "cpu"], default="gpu")
    ap.add_argument("--gpu", help="GPU index, or comma-separated indices on this host to optimize "
                    "a motif's states concurrently (one state per GPU)")
    ap.add_argument("--threads", type=int, default=2)
    ap.add_argument("--nconf", type=int, default=3)
    ap.add_argument("--cpus", help="comma-separated CPU indices; otherwise choose unallocated CPUs")
    ap.add_argument("--thermal-definition", choices=["native", "legacy"], default="native",
                    help="new computations use derivative-checked solution DFT thermochemistry")
    a = ap.parse_args(argv)
    try:
        gpus = [int(v) for v in a.gpu.split(",")] if a.gpu is not None else []
    except ValueError:
        ap.error("--gpu takes integer indices")
    if len(set(gpus)) != len(gpus):
        ap.error("--gpu indices must be distinct (one task per GPU)")
    a.gpu = gpus[0] if gpus else None
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
    workers = max(1, len(gpus))
    if (a.threads * workers > min(len(allowed_cpus), max(1, (os.cpu_count() or 1) // 2))
            or os.getloadavg()[0] + a.threads * workers > (os.cpu_count() or 1)):
        ap.error("insufficient CPU capacity on this shared node")
    taskset = shutil.which("taskset")
    if not taskset:
        ap.error("taskset is required to constrain worker CPU affinity")
    try:
        requested = [int(v) for v in a.cpus.split(",")] if a.cpus else None
        with allocation_guard():
            allocated = cpu_allocation(a.threads * workers, requested)
            node_preflight = host_capacity(a.threads * workers)
            if a.backend == "gpu":
                node_preflight["gpu"] = gpu_suitability(a.gpu)
                node_preflight["gpus"] = [gpu_suitability(g) for g in gpus]
                host = socket.gethostname().split(".")[0]
                available = free_slots([host])
                if any((host, g) not in available for g in gpus):
                    raise RuntimeError("GPU occupied or reserved; no candidate calculation started")
            original_affinity = set(os.sched_getaffinity(0))
            os.sched_setaffinity(0, allocated)
    except (ValueError, RuntimeError, OSError, subprocess.SubprocessError) as exc:
        ap.error(str(exc))
    # Serial stages (and native thermal, which requires exactly --threads CPUs) use the first
    # block; concurrent state workers each get a disjoint block of --threads CPUs.
    cpu_blocks = [allocated[k * a.threads:(k + 1) * a.threads] for k in range(workers)]
    cpu_set = cpu_blocks[0]
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
        manifest = dict(config=str(config), input_sha256=hashlib.sha256(config.read_bytes()).hexdigest(),
                        selected_definitions_sha256=selection_hash(scaffold, groups),
                        ids=ids, host=socket.gethostname(), pid=os.getpid(), python=sys.executable,
                        backend=a.backend, gpu=a.gpu, gpus=gpus, threads=a.threads, cpu_affinity=cpu_set,
                        worker_cpu_affinity=cpu_blocks,
                        node_preflight=node_preflight,
                        xtb=xtb, nconf=a.nconf, thermal_definition=a.thermal_definition,
                        prior_run_manifest=str(attempt_output / "previous_run_manifest.json") if previous else None,
                        attempt_output=str(attempt_output), steps=[])
        mp = output / "run_manifest.json"
        def wrap(argv, gpu, cpus, compute=True):
            argv = [taskset, "--cpu-list", ",".join(map(str, cpus)), *argv]
            return [str(GATE), "run", str(gpu), "--", *argv] if compute and a.backend == "gpu" else argv

        from redox.core.common import read_manifest as declared_rows
        ferrocene_states = tuple(r["state"] for r in declared_rows() if r["id"] == "ferrocene")
        nodes = build_graph(config, groups, attempt_output, a.backend, a.threads, a.nconf,
                            a.thermal_definition, ferrocene_states=ferrocene_states)
        # Review-table building is off the campaign critical path. The parallel driver
        # (more than one GPU) ends after per-record QC; a single-GPU run keeps building
        # tables for existing single-candidate callers. Either way run the full table
        # audit separately with: python -m redox.reporting.tables --out <run>/tables --workers N
        effective_skip_review = a.skip_review or workers > 1
        entries, results = {}, {}
        lockdir = ROOT / "logs/pipeline_locks"

        # Only the coordinator thread (started/finished) mutates manifest and writes it;
        # run() executes in a worker thread and stashes its outcome under its own node
        # name, so no worker ever mutates a dict that atomic_write_json is serializing.
        def started(node):
            log = attempt_output / f"{node.name.replace(':', '_').replace('/', '_')}.log"
            entry = dict(stage=node.label, node=node.name, argv=list(node.argv),
                         log=str(log), started=time.time())
            entries[node.name] = entry
            manifest["steps"].append(entry)
            atomic_write_json(mp, manifest)
            print(f"{node.name}: {log}", flush=True)

        def run(node, slot):
            log = attempt_output / f"{node.name.replace(':', '_').replace('/', '_')}.log"
            if node_done(node):
                results[node.name] = dict(skipped=True, returncode=0, attempts=[])
                return 0
            gpu = gpus[slot] if gpus else a.gpu
            cpus = cpu_blocks[slot]
            argv = wrap(list(node.argv), gpu, cpus, node.compute)
            lock_handle = None
            try:
                if node.lock:
                    # Serialize the one shared ferrocene namespace across any concurrent
                    # runs; the per-gid flocks do not cover it. Blocking, not fail-fast,
                    # so a sibling run waits rather than aborting.
                    lock_handle = (lockdir / f"{node.lock}.lock").open("a+")
                    fcntl.flock(lock_handle, fcntl.LOCK_EX)
                    if node_done(node):
                        results[node.name] = dict(skipped=True, returncode=0, attempts=[])
                        return 0
                proc, attempts = run_stage(argv, log, env)
            finally:
                if lock_handle is not None:
                    fcntl.flock(lock_handle, fcntl.LOCK_UN); lock_handle.close()
            results[node.name] = dict(argv=argv, gpu=gpu, cpu_affinity=cpus,
                                      returncode=proc.returncode, attempts=attempts,
                                      log=attempts[-1]["log"])
            return proc.returncode

        def finished(node, code):
            entries[node.name].update(results.get(node.name, {}), finished=time.time())
            atomic_write_json(mp, manifest)

        result = execute_graph(nodes, run, workers=workers, started=started, finished=finished)
        atomic_write_json(mp, manifest)
        if result["failed"]:
            from redox.quality.readiness import inspect_states
            from redox.core.common import read_manifest
            failed_entry = entries[result["failed"][0]]
            refused = bool(failed_entry.get("attempts")) and failed_entry["attempts"][-1]["reservation_refused"]
            manifest.update(status="resource_unavailable" if refused else "failed",
                            failed_stage=failed_entry["stage"],
                            scientific_accuracy_validated=False,
                            selected_state_qc=inspect_states([r for r in read_manifest() if r["id"] in ids]))
            atomic_write_json(mp, manifest)
            raise SystemExit(f"{failed_entry['node']} failed; see {failed_entry.get('log')}")
        # Written last, after the independent audit passes; QC limits remain in tables.
        from redox.quality.readiness import inspect_states
        from redox.core.common import read_manifest
        states = inspect_states([r for r in read_manifest() if r["id"] in ids])
        if effective_skip_review:
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
            print(f"Verified selected records: {mp}; tables await separate collection "
                  f"(python -m redox.reporting.tables --out {attempt_output / 'tables'} --workers N)")
            return
        # Single-GPU default: build and audit tables as the final, off-DAG step.
        review_log = attempt_output / "review_tables.log"
        review_argv = wrap([sys.executable, "-m", "redox.reporting.tables",
                            "--out", str(attempt_output / "tables")], a.gpu, cpu_set, False)
        entry = dict(stage="review_tables", argv=review_argv, log=str(review_log), started=time.time())
        manifest["steps"].append(entry); atomic_write_json(mp, manifest)
        proc, attempts = run_stage(review_argv, review_log, env)
        entry.update(returncode=proc.returncode, finished=time.time(), attempts=attempts)
        atomic_write_json(mp, manifest)
        if proc.returncode:
            manifest.update(status="failed", failed_stage="review_tables",
                            scientific_accuracy_validated=False, selected_state_qc=states)
            atomic_write_json(mp, manifest)
            raise SystemExit(f"review_tables failed; see {review_log}")
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
