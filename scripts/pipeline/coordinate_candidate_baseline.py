#!/usr/bin/env python
"""Run the candidate's three checked baseline cases and collect a matched cycle.

This is a single-minimum feasibility demonstration alongside broader searches.
Resource and physical pilots are real calculations; no failed input is retried.
"""
import argparse
import fcntl
import json
import os
from pathlib import Path
import signal
import socket
import subprocess
import sys
import time

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.core.protocol import load_record
from redox.qm.native_thermal import load_native_thermal
from redox.workflows import polymer_campaign as api
from redox.workflows.resources import allocation_guard, cpu_allocation, host_capacity

REFERENCE = ROOT / "results/polymer_models/campaigns/2026-10-23_baseline_two_pose_ensemble_v2"
WORKSPACE_DRIVER = ROOT / "scripts/pipeline/workspace_polymer_campaign.py"
XC_WORKSPACE_DRIVER = ROOT / "scripts/pipeline/xc_workspace_campaign.py"
PREPARATION_AFFINITY = tuple(sorted(os.sched_getaffinity(0)))


def preparation_command(command):
    """Let the preparation CLI allocate separately from the pinned observer.

    The CLI's normal allocation guard then selects/pins its own free CPU. A
    child confined to its observer parent's one occupied CPU cannot do that.
    Scientific worker affinities and thread limits are unchanged.
    """
    return ["taskset", "--cpu-list", ",".join(map(str, PREPARATION_AFFINITY)), *command]


def same_model(parent, replacement):
    for field in ("quantum_settings", "thermochemistry", "solvent"):
        if parent[field] != replacement[field]:
            raise ValueError("Execution recovery changed the chemical/numerical model")
    if parent["tasks"] != replacement["tasks"]:
        raise ValueError("Execution recovery changed the exact task identities")


def resolved_declaration(base, declaration):
    path = base / "baseline_recoveries.json"
    if not path.exists():
        return declaration
    registry = json.loads(path.read_text())
    cases = []
    for case in declaration["cases"]:
        replacement = registry.get("overrides", {}).get(case["key"])
        if not replacement:
            cases.append(case)
            continue
        parent, target = Path(case["campaign"]), Path(replacement["campaign"])
        original_plan, new_plan = api.read_plan(parent), api.read_plan(target)
        same_model(original_plan, new_plan)
        api.verify_sources(new_plan["producer_sha256"])
        api.verify_sources(new_plan["input_sha256"])
        cases.append(dict(case, campaign=str(target), recovered_from=str(parent),
            recovery_scope=new_plan["recovery_scope"], execution_hosts=replacement.get("hosts")))
    return dict(declaration, cases=cases)


def reconcile_external_energy(base):
    """Register native work only after an external electronic parent is checked.

    The source calculation and any existing target are preserved. Creating a
    thermal workflow does not create a native result or repeat an optimization.
    """
    registry_path = base / "baseline_recoveries.json"
    if not registry_path.exists():
        return []
    events = []
    with registry_path.with_suffix(".lock").open("a+") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        registry = json.loads(registry_path.read_text())
        external = registry.get("external_continuations", {})
        for key, item in list(external.items()):
            if not item.get("native_target"):
                continue
            source = Path(item["campaign"])
            source = source if source.is_absolute() else ROOT / source
            workflow = source / "workspace_energy_status.json"
            if not workflow.exists() or json.loads(workflow.read_text()).get("status") != "stationary_electronic_case_checked_native_G_pending":
                continue
            target = Path(item["native_target"])
            target = target if target.is_absolute() else ROOT / target
            controls = Path(item["native_controls"])
            controls = controls if controls.is_absolute() else ROOT / controls
            execution_driver = Path(item.get("native_execution_driver", str(WORKSPACE_DRIVER))).resolve()
            if execution_driver not in (WORKSPACE_DRIVER, XC_WORKSPACE_DRIVER):
                raise ValueError("Unregistered native handoff execution driver")
            module = api.load_module("candidate_native_handoff", WORKSPACE_DRIVER)
            if not target.exists():
                if execution_driver == WORKSPACE_DRIVER:
                    module.prepare_native(source.resolve(), target.resolve(), controls.resolve(),
                        item["native_workspace_GiB"], item["native_minimum_GPU_MiB"])
                else:
                    subprocess.run(preparation_command([sys.executable, str(execution_driver), "prepare-native", "--source", str(source.resolve()),
                        "--out", str(target.resolve()), "--controls", str(controls.resolve()),
                        "--workspace-GiB", str(item["native_workspace_GiB"]),
                        "--minimum-GPU-MiB", str(item["native_minimum_GPU_MiB"])]), check=True)
            original, new = api.read_plan(source), api.read_plan(target)
            same_model(original, new)
            if (new.get("source_campaign") != str(source.resolve())
                    or not new.get("inherited_energy_only")
                    or new.get("execution_driver") != str(execution_driver)):
                raise ValueError("Existing native target belongs to a different handoff")
            api.verify_sources(new["producer_sha256"])
            api.verify_sources(new["input_sha256"])
            records, invalid = api.admitted_records(target, "energy")
            if invalid or len(records) != 1 or next(iter(records.values()))["task"]["model_key"] != key:
                raise ValueError("Native handoff lacks its exact checked electronic parent")
            history = base / "baseline_observer_history"
            history.mkdir(exist_ok=True)
            previous = history / f"registry_before_native_handoff_{time.time_ns()}.json"
            previous.write_bytes(registry_path.read_bytes())
            registry.setdefault("overrides", {})[key] = dict(campaign=str(target.resolve()),
                reason="Exact checked electronic parent reused; separate native-thermal resource and physical qualification pending")
            if item.get("native_hosts"):
                registry["overrides"][key]["hosts"] = item["native_hosts"]
            external.pop(key)
            registry.update(previous_registry=str(previous), updated=time.time())
            atomic_write_json(registry_path, registry)
            events.append(dict(key=key, source=str(source.resolve()), target=str(target.resolve()),
                registered=time.time(), repeated_optimization=False, native_G_computed_by_handoff=False))
    if events:
        with (base / "native_handoff_events.jsonl").open("a") as stream:
            for event in events:
                stream.write(json.dumps(event) + "\n")
    return events


def collect(base, declaration, driver):
    if len(declaration["cases"]) != 3 or {tuple(row["occupied_sites"]) for row in declaration["cases"]} != {(2,), (3,), (2, 3)}:
        raise ValueError("The candidate baseline needs exactly the three declared adjacent-site components")
    rows = []
    for case in declaration["cases"]:
        report = driver.collect(Path(case["campaign"]))
        if len(report["records"]) != 1:
            raise ValueError("The candidate baseline requires one explicitly declared start per occupied graph")
        row = report["records"][0]
        rows.append(dict(case, E_eV=row["E_eV"], G_eV=row["G_eV"], checked_energy=report["checked_energy"],
            checked_native_G=report["checked_native_G"], conformer_convergence_proven=False))
    references, invalid = api.admitted_records(REFERENCE, "thermal")
    if invalid:
        raise ValueError("Shared empty reference contains invalid receipts")
    matches = [row for row in references.values() if row["task"]["smiles"] == declaration["shared_empty_reference"]["graph_smiles"]]
    if len(matches) != 2:
        raise ValueError("Both checked empty-backbone references must remain available")
    reference = min(matches, key=lambda row: row["E_smd_eV"])
    state = Path(reference["state_dir"])
    native = load_native_thermal(state, 0, 1)
    energy = load_record(state, 0, 1)
    if native is None or energy is None or abs(energy["e_smd_eV"] + native["g_thermal_total_eV"] - reference["G_eV"]) > 1e-8:
        raise ValueError("Shared empty-component G failed ordinary native admission")
    terms = json.loads((base / "reaction_terms.json").read_text())
    api.verify_sources(terms["source_sha256"])
    lookup = {tuple(row["occupied_sites"]): row for row in rows}
    cycle = None
    if all(row["G_eV"] is not None for row in rows):
        empty, first, second, both = reference["G_eV"], lookup[(2,)]["G_eV"], lookup[(3,)]["G_eV"], lookup[(2, 3)]["G_eV"]
        penalty = (both + empty - first - second) * 96.48533212331001
        models = {}
        for name, constant in terms["constant_eV"].items():
            if constant is None:
                models[name] = None
                continue
            models[name] = dict(graft_3_on_empty_kJ_mol=(second - empty + constant) * 96.48533212331001,
                graft_3_given_2_kJ_mol=(both - first + constant) * 96.48533212331001,
                graft_2_on_empty_kJ_mol=(first - empty + constant) * 96.48533212331001,
                graft_2_given_3_kJ_mol=(both - second + constant) * 96.48533212331001)
        cycle = dict(deltaDeltaG_kJ_mol=penalty, DBU_models=models,
            interpretation="Neutral-chain conditional site nonadditivity; not isolated motif binding, electronic coupling or an ET rate")
    candidate = json.loads((base / "candidate_snapshot.json").read_text())
    routes = json.loads((base / "route_snapshot.json").read_text())
    current = dict(status="single_minimum_baseline_components_complete" if cycle is not None else "candidate_baseline_components_pending",
        candidate_id=candidate["id"], exact_capped_product_key=json.loads((base / "plan.json").read_text())["exact_capped_product_key"],
        cases=rows, shared_empty_component=dict(state_dir=str(state), G_eV=reference["G_eV"],
            selection="Lowest electronic-energy member of the two checked exact empty-graph states; no occupied control energy transferred"),
        conditional_polymer_cycle=cycle, source_capped_redox=dict(E1_vs_Fc_V=candidate.get("E1_vs_Fc_V"), E2_vs_Fc_V=candidate.get("E2_vs_Fc_V"),
            scope="Exact capped graft only; not an assigned polymer potential"),
        prospective_route_status=routes["gate"]["route_status"], missing_route_evidence=routes["gate"]["missing_route_evidence"],
        accessible_material_capacity_established=False, synthesis_admission=False,
        broader_sampling_complete=False, updated=time.time())
    atomic_write_json(base / "baseline_current.json", current)
    return current


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--campaign", type=Path, required=True)
    parser.add_argument("--hosts", default="lambda13")
    parser.add_argument("--max-hours", type=float, default=48)
    parser.add_argument("--resume-observer", action="store_true")
    args = parser.parse_args()
    if not 0 < args.max_hours <= 48:
        parser.error("Use a finite duration up to 48 hours")
    base = args.campaign.resolve()
    declaration = json.loads((base / "baseline_plan.json").read_text())
    api.verify_sources(declaration["source_sha256"])
    shared = Path(declaration["resource_root"])
    driver_path = ROOT / "scripts/pipeline/refine_polymer_conformers.py"
    orchestration_path = ROOT / "scripts/pipeline/coordinate_crest_refinement.py"
    driver = api.load_module("candidate_baseline_QM", driver_path)
    orchestration = api.load_module("candidate_baseline_launch_control", orchestration_path)
    sources = {str(p): api.digest(p) for p in (Path(__file__).resolve(), driver_path, orchestration_path,
        WORKSPACE_DRIVER, XC_WORKSPACE_DRIVER, base / "baseline_plan.json")}
    with allocation_guard():
        host_capacity(1)
        os.sched_setaffinity(0, cpu_allocation(1))
    with (base / "baseline_coordinator.lock").open("a+") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        path = base / "baseline_coordinator.json"
        if path.exists():
            old = json.loads(path.read_text())
            ownership = api.load_module("candidate_observer_ownership", ROOT / "scripts/pipeline/run_polymer_spacing.py")
            if not args.resume_observer or ownership.owner_alive(old):
                raise RuntimeError("Existing candidate observer is live or explicit resume was not requested")
            history = base / "baseline_observer_history"
            history.mkdir(exist_ok=True)
            (history / (str(time.time_ns()) + ".json")).write_bytes(path.read_bytes())
        deadline = old["deadline"] if path.exists() and args.resume_observer else time.time() + args.max_hours * 3600
        owner = dict(host=socket.gethostname(), pid=os.getpid(), process_start=Path("/proc/self/stat").read_text().rsplit(")", 1)[1].split()[19],
            started=time.time(), deadline=deadline, source_sha256=sources,
            CPU_affinity=sorted(os.sched_getaffinity(0)), status="observing_submitted_energy_pilots")
        def stop(signum, frame):
            raise KeyboardInterrupt(f"signal {signum}; scientific workers preserved")
        signal.signal(signal.SIGTERM, stop)
        atomic_write_json(path, owner)
        try:
            while time.time() < owner["deadline"]:
                api.verify_sources(sources)
                reconcile_external_energy(base)
                active = resolved_declaration(base, declaration)
                for case in active["cases"]:
                    campaign = Path(case["campaign"])
                    report = driver.collect(campaign)
                    for stage in ("energy", "thermal"):
                        # A claim may be live or its ownership may be unresolved.
                        # Do not start extra workers just to discover the same
                        # claimed single task. Stale-owner recovery is separate.
                        if any((campaign / stage / "claims").iterdir()):
                            continue
                        if any((campaign / stage / "failures").glob("*.json")):
                            continue
                        if stage == "thermal" and (report["checked_energy"] != 1 or report["checked_native_G"] == 1):
                            continue
                        if stage == "energy" and report["checked_energy"] == 1:
                            continue
                        pilot_path = campaign / stage / "pilot/status.json"
                        pilot = json.loads(pilot_path.read_text()) if pilot_path.exists() else {}
                        if pilot.get("status") == "failed_preserved":
                            continue
                        original_driver = orchestration.DRIVER
                        orchestration.DRIVER = Path(api.read_plan(campaign).get("execution_driver", original_driver))
                        try:
                            if pilot.get("status") in ("memory_pilot_passed", "thermal_pilot_passed"):
                                orchestration.launch(campaign, stage, False, case.get("execution_hosts") or args.hosts.split(","), shared)
                            elif not any((campaign / stage / "launches").glob("*_pilot.json")):
                                orchestration.launch(campaign, stage, True, case.get("execution_hosts") or args.hosts.split(","), shared)
                        finally:
                            orchestration.DRIVER = original_driver
                current = collect(base, active, driver)
                owner.update(updated=time.time(), status=current["status"], checked_energy=sum(row["checked_energy"] for row in current["cases"]),
                    checked_native_G=sum(row["checked_native_G"] for row in current["cases"]))
                atomic_write_json(path, owner)
                if current["conditional_polymer_cycle"] is not None:
                    return
                time.sleep(30)
            owner["status"] = "finite_candidate_baseline_window_ended_outputs_preserved"
        except BaseException as error:
            owner.update(status="candidate_baseline_observation_needs_review_outputs_preserved", error=repr(error))
            raise
        finally:
            owner["finished"] = time.time()
            atomic_write_json(path, owner)


if __name__ == "__main__":
    main()
