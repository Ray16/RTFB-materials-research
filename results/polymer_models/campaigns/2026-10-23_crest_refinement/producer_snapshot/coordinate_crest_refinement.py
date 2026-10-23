#!/usr/bin/env python
"""Bounded final-CREST -> DFT -> native-G continuation and matched cycles.

No source search is restarted and no failed state is retried by this controller.
Energy cases can overlap; one native case at a time shares its measured RAM
budget. The separate baseline native panel must release its live workers first.
"""
import argparse
import fcntl
import json
import os
from pathlib import Path
import shlex
import signal
import socket
import subprocess
import sys
import time

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api
from redox.workflows.resources import allocation_guard, cpu_allocation, host_capacity

SCRIPT = Path(__file__).resolve()
DRIVER = ROOT / "scripts/pipeline/refine_polymer_conformers.py"
PREPARER = ROOT / "scripts/pipeline/prepare_crest_refinement.py"
BASELINE = ROOT / "results/polymer_models/campaigns/2026-10-23_baseline_two_pose_ensemble_v2"


def observe_workers(root):
    """One ownership RPC per host; unreachable hosts are unknown, never dead."""
    paths = list(root.glob("*/thermal/workers/*.json")) + list(root.glob("*/energy/workers/*.json"))
    paths += list((BASELINE / "thermal/workers").glob("*.json"))
    by_host = {}
    for path in paths:
        row = json.loads(path.read_text())
        by_host.setdefault(row["host"], []).append(dict(path=str(path), pid=row["pid"], start=row["process_start"],
            stage=row["stage"], baseline=path.is_relative_to(BASELINE), case=path.parents[2].name))
    observed = []
    code = '''import json,sys
from pathlib import Path
rows=json.load(sys.stdin)
for row in rows:
 p=Path('/proc')/str(row['pid'])/'stat'
 row['live']=p.exists() and p.read_text().rsplit(')',1)[1].split()[19]==row['start']
 row['unknown']=False
print(json.dumps(rows))'''
    for host, rows in by_host.items():
        try:
            result = subprocess.run(["ssh", "-o", "BatchMode=yes", "-o", "ConnectTimeout=8", host, shlex.join([sys.executable, "-c", code])],
                input=json.dumps(rows), capture_output=True, text=True, timeout=20, check=True)
            observed.extend(json.loads(result.stdout))
        except (OSError, ValueError, subprocess.SubprocessError):
            observed.extend(dict(row, live=None, unknown=True) for row in rows)
    return observed


def launch(campaign, stage, pilot, hosts, root):
    for host in hosts:
        preflight = root / (host + "_preflight.json")
        data = json.loads(preflight.read_text())
        if data.get("ok") is not True or data.get("python") != sys.executable:
            raise ValueError("Required node/interpreter preflight is missing")
        target = campaign / (host + "_preflight.json")
        if not target.exists():
            target.write_bytes(preflight.read_bytes())
        argv = ["env", "PYTHONPATH=" + str(ROOT / "src"),
            *[key + "=1" for key in ("OMP_NUM_THREADS", "OMP_THREAD_LIMIT", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS", "XTB_THREADS")],
            sys.executable, str(DRIVER), "launch", "--campaign", str(campaign), "--stage", stage]
        if pilot:
            argv.append("--pilot")
        with (root / "launch_observations.log").open("a") as log:
            result = subprocess.run(["ssh", "-n", "-o", "BatchMode=yes", "-o", "ConnectTimeout=8", host, shlex.join(argv)], stdout=log, stderr=subprocess.STDOUT, timeout=55)
        if result.returncode:
            raise RuntimeError("Launch observation needs review; no calculation failure inferred")
        # Pilots are serial dependencies, not one pilot per host.
        if pilot and any((campaign / stage / "launches").glob("*_pilot.json")):
            break


def matched_cycles(request, reports):
    """Join only matching constitutional cases within one backbone configuration."""
    from redox.properties.conformer_ensemble import neighbor_cycle

    lookup = {}
    for entry in request["cases"]:
        case = json.loads(Path(entry["case_path"]).read_text())
        lookup[(entry["configuration"], tuple(case["occupied_sites"]))] = entry["key"]
        report = reports.get(entry["key"])
        if report and report.get("primary_G_eV") is not None:
            if not report.get("partition_protocol") or report.get("temperature_K") != request["temperature_K"] or report.get("solvent") != request["solvent"]:
                raise ValueError("Cycle component lacks its declared thermodynamic context")
            for field in ("smiles", "backbone_stereo_fingerprint", "occupied_sites", "reactive_sites", "repeat_units"):
                if report["case"].get(field) != case.get(field):
                    raise ValueError("Cycle component belongs to another exact occupancy/stereo graph")
    cycles = []
    for configuration in sorted({key[0] for key in lookup}):
        for sites in ((2, 3), (2, 4)):
            keys = [lookup[(configuration, occupied)] for occupied in ((), (sites[0],), (sites[1],), sites)]
            values = [reports.get(key, {}).get("primary_G_eV") for key in keys]
            protocols = [reports[key].get("partition_protocol") for key in keys if key in reports and reports[key].get("primary_G_eV") is not None]
            if protocols and any(value != protocols[0] for value in protocols):
                raise ValueError("Matched cycle combines different thermal partitions")
            cycles.append(dict(configuration=configuration, sites=list(sites), case_keys=keys,
                deltaDeltaG_kJ_mol=neighbor_cycle(*values) if all(value is not None for value in values) else None,
                scope="Conditional finite-set site nonadditivity; no experimental tacticity weights or pure motif-binding interpretation"))
    return cycles


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--hosts", required=True, help="Comma-separated preflighted suitable hosts")
    parser.add_argument("--max-hours", type=float, default=48)
    args = parser.parse_args()
    if not 0 < args.max_hours <= 48:
        parser.error("Use a finite controller window up to 48 hours")
    hosts = args.hosts.split(",")
    if any(not host or host.startswith("-") for host in hosts):
        parser.error("Invalid node list")
    root = args.root.resolve()
    source = {str(p): api.digest(p) for p in (SCRIPT, DRIVER, PREPARER, root / "request.json")}
    with allocation_guard():
        host_capacity(1)
        os.sched_setaffinity(0, cpu_allocation(1))
    with (root / "coordinator.lock").open("a+") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        path = root / "coordinator.json"
        if path.exists():
            raise RuntimeError("Prior controller ownership must be reviewed; no blind restart")
        owner = dict(host=socket.gethostname(), pid=os.getpid(), process_start=Path("/proc/self/stat").read_text().rsplit(")", 1)[1].split()[19],
            started=time.time(), deadline=time.time() + args.max_hours * 3600, source_sha256=source, CPU_affinity=sorted(os.sched_getaffinity(0)), status="starting")
        def stop(signum, frame):
            raise KeyboardInterrupt(f"signal {signum}; submitted quantum workers preserved")
        signal.signal(signal.SIGTERM, stop)
        atomic_write_json(path, owner)
        previous = {}
        request = json.loads((root / "request.json").read_text())
        preparer = api.load_module("final_CREST_stager", PREPARER)
        driver = api.load_module("native_refinement_collector", DRIVER)
        try:
            while time.time() < owner["deadline"]:
                api.verify_sources(source)
                preparation = preparer.stage(root)
                live = observe_workers(root)
                native_busy = [row for row in live if row["stage"] == "thermal" and (row["live"] or row["unknown"])]
                active_native_case = next((row["case"] for row in native_busy if not row["baseline"]), None)
                baseline_busy = any(row["baseline"] for row in native_busy)
                baseline_controller = BASELINE / "coordinator.json"
                baseline_state = json.loads(baseline_controller.read_text()) if baseline_controller.exists() else {}
                # A completed pilot can leave a brief gap before the baseline
                # controller launches its queue. Do not mistake that gap for a
                # release of the panel's future host-memory budget.
                baseline_terminal = baseline_state.get("status") in {
                    "finite_two_pose_partitions_complete", "pilot_failed_outputs_preserved",
                    "finite_native_panel_terminal_with_explicit_failures", "finite_window_ended_outputs_preserved",
                    "observation_or_workflow_needs_review_outputs_preserved",
                } and "finished" in baseline_state
                native_uncertain = any(row["unknown"] for row in native_busy)
                reports = {}
                for entry in request["cases"]:
                    campaign = root / entry["key"]
                    if not (campaign / "plan.json").exists():
                        continue
                    inventory = tuple((str(p), api.digest(p)) for phase in ("energy", "thermal") for p in sorted((campaign / phase / "done").glob("*.json")))
                    current_path = campaign / "current.json"
                    if inventory != previous.get(entry["key"]) or not current_path.exists():
                        report = driver.collect(campaign)
                        previous[entry["key"]] = inventory
                    else:
                        report = json.loads(current_path.read_text())
                    reports[entry["key"]] = report
                    for stage in ("energy", "thermal"):
                        if any((campaign / stage / "failures").glob("*.json")):
                            continue
                        pilot_path = campaign / stage / "pilot/status.json"
                        pilot = json.loads(pilot_path.read_text()) if pilot_path.exists() else {}
                        if pilot.get("status") == "failed_preserved":
                            continue
                        if stage == "thermal":
                            if report["checked_energy"] != report["expected_starts"] or baseline_busy or not baseline_terminal or native_uncertain:
                                continue
                            if active_native_case not in (None, entry["key"]):
                                continue
                            if report["checked_native_G"] == report["expected_starts"]:
                                continue
                            active_native_case = entry["key"]
                        elif report["checked_energy"] == report["expected_starts"]:
                            continue
                        if pilot.get("status") in ("memory_pilot_passed", "thermal_pilot_passed"):
                            launch(campaign, stage, False, hosts, root)
                        elif not any((campaign / stage / "launches").glob("*_pilot.json")):
                            launch(campaign, stage, True, hosts, root)
                cycles = matched_cycles(request, reports)
                current = dict(preparation=preparation, cases=reports, cycles=cycles,
                    observed_worker_ownership=live, updated=time.time(),
                    finite_partitions_complete=sum(report["primary_G_eV"] is not None for report in reports.values()),
                    expected_cases=len(request["cases"]), global_sampling_convergence_proven=False,
                    synthesis_candidates_admitted=0)
                atomic_write_json(root / "current.json", current)
                owner.update(status="source_sampling_DFT_native_followup_in_progress", updated=time.time(),
                    frozen_graphs=preparation["frozen_graphs"], checked_energy=sum(report["checked_energy"] for report in reports.values()),
                    checked_native_G=sum(report["checked_native_G"] for report in reports.values()))
                atomic_write_json(path, owner)
                if current["finite_partitions_complete"] == current["expected_cases"]:
                    owner["status"] = "finite_initial_partitions_complete_sampling_diagnostics_require_review"
                    return
                time.sleep(30)
            owner["status"] = "finite_controller_window_ended_workers_and_outputs_preserved"
        except BaseException as error:
            owner.update(status="workflow_or_observation_needs_review_outputs_preserved", error=repr(error))
            raise
        finally:
            owner["finished"] = time.time()
            atomic_write_json(path, owner)


if __name__ == "__main__":
    main()
