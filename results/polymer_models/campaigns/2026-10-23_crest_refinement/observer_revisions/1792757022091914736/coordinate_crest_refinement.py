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
LAUNCH_REMOTE = '''import json,os,subprocess,sys,time
from pathlib import Path
c=json.load(sys.stdin);handle=Path(c['handle']);outcome=Path(c['outcome'])
if handle.exists():
 old=json.loads(handle.read_text());p=Path('/proc')/str(old['pid'])/'stat'
 live=False
 try:
  fields=p.read_text().rsplit(')',1)[1].split();live=fields[19]==old['process_start'] and fields[0]!='Z'
 except FileNotFoundError:pass
 if live:
  print(json.dumps(dict(status='same_launch_observer_live',handle=old)));sys.exit(0)
 report=json.loads(outcome.read_text()) if outcome.exists() else {}
 if report.get('pid')!=old['pid'] or report.get('process_start')!=old['process_start'] or 'finished' not in report:
  print(json.dumps(dict(status='ended_launch_outcome_unconfirmed_preserve',handle=old)));sys.exit(0)
 if report.get('status')!='launch_observation_completed_outputs_preserved':
  print(json.dumps(dict(status='launch_outcome_needs_review_preserve',handle=old)));sys.exit(0)
 with (handle.parent/'launch_observer_history.jsonl').open('a') as stream:stream.write(json.dumps(dict(handle=old,outcome=report))+'\\n')
env=dict(os.environ,PYTHONPATH=c['root']+'/src')
for key in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS','XTB_THREADS'):env[key]='1'
with open(c['log'],'a') as log:
 child=subprocess.Popen(c['argv'],cwd=c['root'],env=env,stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
record=dict(host=os.uname().nodename,pid=child.pid,process_start=Path(f'/proc/{child.pid}/stat').read_text().rsplit(')',1)[1].split()[19],argv=c['argv'],started=time.time())
temporary=handle.with_suffix('.tmp');temporary.write_text(json.dumps(record,indent=2)+'\\n');temporary.replace(handle)
print(json.dumps(dict(status='launch_observer_started_scientific_starts_unconfirmed',handle=record)))'''


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
        argv = [sys.executable, str(DRIVER), "launch", "--campaign", str(campaign), "--stage", stage]
        if pilot:
            argv.append("--pilot")
        role = "pilot" if pilot else "queue"
        payload = dict(root=str(ROOT), argv=argv,
            handle=str(campaign / stage / f"launcher_handle_{host}_{role}.json"),
            outcome=str(campaign / stage / f"launch_observer_{host}_{role}.json"),
            log=str(campaign / stage / f"launch_observer_{host}_{role}.log"))
        result = subprocess.run(["ssh", "-o", "BatchMode=yes", "-o", "ConnectTimeout=8", host,
            shlex.join([sys.executable, "-c", LAUNCH_REMOTE])], input=json.dumps(payload), capture_output=True, text=True, timeout=25, check=True)
        observation = json.loads(result.stdout)
        with (root / "launch_observations.log").open("a") as log:
            log.write(json.dumps(dict(campaign=str(campaign), stage=stage, role=role, observation=observation, checked=time.time())) + "\n")
        # Pilots are serial dependencies, not one pilot per host.
        if pilot:
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
    parser.add_argument("--resume-observer", action="store_true", help="Resume only after proving the previous control process has ended")
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
            previous_owner = json.loads(path.read_text())
            owner_check = api.load_module("cohort_observer_ownership", ROOT / "scripts/pipeline/run_polymer_spacing.py")
            if not args.resume_observer or owner_check.owner_alive(previous_owner):
                raise RuntimeError("Prior controller is still live or explicit observer resume was not requested")
            history = root / "coordinator_history"
            history.mkdir(exist_ok=True)
            (history / (str(time.time_ns()) + ".json")).write_bytes(path.read_bytes())
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
                baseline_current = json.loads((BASELINE / "current.json").read_text())
                atomic_write_json(BASELINE / "continuation_status.json", dict(
                    status="native_ownership_temporarily_unverified" if any(row["unknown"] for row in native_busy if row["baseline"]) else "observing_existing_native_workers" if baseline_busy else "native_panel_complete" if baseline_current["checked_native_G"] == baseline_current["expected_DFT_poses"] else "no_live_native_handle_requires_review",
                    checked_native_G=baseline_current["checked_native_G"], expected_native_G=baseline_current["expected_DFT_poses"],
                    observed_workers=[row for row in live if row["baseline"]], observer_owner=str(path),
                    scientific_jobs_restarted=0, updated=time.time()))
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
                    baseline_native_G=baseline_current["checked_native_G"],
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
