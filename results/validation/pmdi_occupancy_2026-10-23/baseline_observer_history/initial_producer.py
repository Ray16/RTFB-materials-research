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
import time

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.core.protocol import load_record
from redox.qm.native_thermal import load_native_thermal
from redox.workflows import polymer_campaign as api
from redox.workflows.resources import allocation_guard, cpu_allocation, host_capacity

REFERENCE = ROOT / "results/polymer_models/campaigns/2026-10-23_baseline_two_pose_ensemble_v2"


def collect(base, declaration, driver):
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
        candidate_id=candidate["id"], exact_capped_product_key=declaration.get("exact_capped_product_key"),
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
    sources = {str(p): api.digest(p) for p in (Path(__file__).resolve(), driver_path, orchestration_path, base / "baseline_plan.json")}
    with allocation_guard():
        host_capacity(1)
        os.sched_setaffinity(0, cpu_allocation(1))
    with (base / "baseline_coordinator.lock").open("a+") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        path = base / "baseline_coordinator.json"
        if path.exists():
            raise RuntimeError("Existing candidate observer ownership requires review")
        owner = dict(host=socket.gethostname(), pid=os.getpid(), process_start=Path("/proc/self/stat").read_text().rsplit(")", 1)[1].split()[19],
            started=time.time(), deadline=time.time() + args.max_hours * 3600, source_sha256=sources,
            CPU_affinity=sorted(os.sched_getaffinity(0)), status="observing_submitted_energy_pilots")
        def stop(signum, frame):
            raise KeyboardInterrupt(f"signal {signum}; scientific workers preserved")
        signal.signal(signal.SIGTERM, stop)
        atomic_write_json(path, owner)
        try:
            while time.time() < owner["deadline"]:
                api.verify_sources(sources)
                for case in declaration["cases"]:
                    campaign = Path(case["campaign"])
                    report = driver.collect(campaign)
                    for stage in ("energy", "thermal"):
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
                        if pilot.get("status") in ("memory_pilot_passed", "thermal_pilot_passed"):
                            orchestration.launch(campaign, stage, False, args.hosts.split(","), shared)
                        elif not any((campaign / stage / "launches").glob("*_pilot.json")):
                            orchestration.launch(campaign, stage, True, args.hosts.split(","), shared)
                current = collect(base, declaration, driver)
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
