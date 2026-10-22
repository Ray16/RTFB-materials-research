#!/usr/bin/env python
"""Finish the neutral spacing cohort with qualified finer CDS derivatives.

Recover the failed large pilot from its saved Hessian, then use the same
qualified derivative setting for remaining states. Preserve prior failures,
borrow checked energies and native results, and retain all physical gates.
"""

import argparse
import json
import os
from pathlib import Path
import shutil
import socket
import time

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api
from redox.workflows.checkpoints import verify_state_checkpoint

SCRIPT = Path(__file__).resolve()
CAMPAIGN = ROOT / "results/polymer_models/campaigns/2026-10-21_spacing_cds_recovery"
PARENT = ROOT / "results/polymer_models/campaigns/2026-10-21_spacing_bounded_pcm"
KEY = "occ_2_4__pose1"


def prepare():
    force_review = (
        ROOT / "results/validation/spacing_CDS_forces_frozen_2026-10-21/summary.json"
    )
    force = json.loads(force_review.read_text())
    if force["status"] != "energy_gradient_diagnostics_passed" or any(
        row["maximum_legacy_disagreement"] > 1e-6
        or row["maximum_step_disagreement"] > 1e-6
        for row in force["cases"]
    ):
        raise ValueError(
            "Selected spacing geometries require conservative force review"
        )
    parent = api.read_plan(PARENT)
    api.verify_sources(parent["producer_sha256"])
    energy, invalid = api.admitted_records(PARENT, "energy")
    repaired = (
        ROOT
        / "results/polymer_models/campaigns/2026-10-21_spacing_singlepoint_repair"
        / "energy/done/occ_none__pose1.json"
    )
    replacement = verify_state_checkpoint(repaired)
    energy[replacement["task"]["key"]] = replacement
    if invalid or len(energy) != 12:
        raise ValueError("All twelve checked electronic conformers required")
    CAMPAIGN.mkdir(exist_ok=False)
    for stage in ("energy", "thermal"):
        for name in (
            "done",
            "tasks",
            "claims",
            "workers",
            "failures",
            "pilot",
            "launches",
        ):
            (CAMPAIGN / stage / name).mkdir(parents=True)
    for key in energy:
        source = (
            repaired
            if key == "occ_none__pose1"
            else PARENT / "energy/done" / (key + ".json")
        )
        shutil.copyfile(source, CAMPAIGN / "energy/done" / source.name)
    selected = {}
    for row in energy.values():
        model = row["task"]["model_key"]
        if model not in selected or row["E_smd_eV"] < selected[model]["E_smd_eV"]:
            selected[model] = row
    for row in force["cases"]:
        if (
            api.digest(Path(selected[row["id"]]["state_dir"]) / "opt.xyz")
            != row["geometry_sha256"]
        ):
            raise ValueError("Force audit is for a different selected geometry")
    previous_G = (
        ROOT
        / "results/polymer_models/campaigns/2026-10-21_spacing_extra_thermal"
        / "thermal/done/occ_none__pose1.json"
    )
    verify_state_checkpoint(previous_G)
    shutil.copyfile(previous_G, CAMPAIGN / "thermal/done" / previous_G.name)
    for name in ("cases.json", "reaction_terms.json", "replay_qualification.json"):
        shutil.copyfile(PARENT / name, CAMPAIGN / name)
    shutil.copytree(PARENT / "controls", CAMPAIGN / "controls")
    shutil.copytree(PARENT / "live_inputs", CAMPAIGN / "live_inputs")
    shutil.copyfile(
        PARENT / "energy/pilot/status.json", CAMPAIGN / "energy/pilot/status.json"
    )
    raw_source = ROOT / f"calcs/polymer_models/{PARENT.name}/thermal/{KEY}/rest"
    original = next(raw_source.glob("native_thermal/*/hessian_cartesian.npy"))
    components = ROOT / "results/validation/polymer_hessian_components_2026-10-21"
    replay = (
        ROOT
        / "results/validation/hessian_reconstruction_replay_2026-10-21/summary.json"
    )
    if (
        json.loads(replay.read_text())["status"]
        != "two_full_recomputation_replays_passed"
    ):
        raise ValueError("Full reconstruction control replays have not passed")
    spec = dict(
        source_hessian=str(original),
        coarse_CDS=str(components / "CDS_hessian_0.0001.npy"),
        fine_CDS=str(components / "CDS_hessian_1e-05.npy"),
        refined_hessian=str(components / "diagnostic_replacement_1e-05.npy"),
        reconstruction_replay=str(replay),
        source_geometry_sha256=api.digest(raw_source / "opt.xyz"),
        memory_replay_evidence_sha256={
            str(p): api.digest(p) for p in (PARENT / "controls").glob("*/status.json")
        },
    )
    files = [
        original,
        components / "CDS_hessian_0.0001.npy",
        components / "CDS_hessian_1e-05.npy",
        components / "diagnostic_replacement_1e-05.npy",
        replay,
        force_review,
        raw_source / "opt.xyz",
        raw_source / "solvent_memory_execution.json",
        original.parent / "producer_source.py",
    ]
    spec["input_sha256"] = {str(p): api.digest(p) for p in files}
    atomic_write_json(CAMPAIGN / "cached_hessian.json", spec)
    plan = dict(
        parent,
        raw_root=str(ROOT / "calcs/polymer_models" / CAMPAIGN.name),
        property_scope=(
            "Same neutral spacing panel; qualified CDS refinement; "
            "no changed physical thresholds"
        ),
    )
    new_sources = [
        SCRIPT,
        ROOT / "src/redox/qm/hessian_resume.py",
        ROOT / "src/redox/qm/cds_finite.py",
        ROOT / "scripts/validation/diagnose_saved_hessian.py",
        ROOT / "scripts/validation/replay_hessian_reconstruction.py",
    ]
    plan["producer_sha256"] = dict(
        parent["producer_sha256"], **{str(p): api.digest(p) for p in new_sources}
    )
    inputs = [
        PARENT / "plan.json",
        repaired,
        previous_G,
        CAMPAIGN / "cached_hessian.json",
        replay,
    ]
    inputs += list((CAMPAIGN / "energy/done").glob("*.json"))
    plan["input_sha256"] = dict(
        parent["input_sha256"], **{str(p): api.digest(p) for p in inputs}
    )
    atomic_write_json(CAMPAIGN / "plan.json", plan)
    api.collect(CAMPAIGN)


def seed_electronic_record(gid, state_name, state, q, m, **kwargs):
    from redox.qm import sp
    from redox.core.protocol import load_record

    selected, _ = api.admitted_records(CAMPAIGN, "energy")
    source = Path(selected[gid]["state_dir"])
    if api.digest(source / "opt.xyz") != api.digest(Path(state) / "opt.xyz"):
        raise ValueError("Electronic cache belongs to a different geometry")
    target = Path(state) / "sp"
    target.mkdir(exist_ok=True)
    for path in (source / "sp").glob("*.json"):
        destination = target / path.name
        if not destination.exists():
            shutil.copyfile(path, destination)
    admitted = load_record(state, q, m)
    if admitted is None:
        raise ValueError("Copied active SP does not pass ordinary provenance admission")
    if abs(admitted["e_smd_eV"] - selected[gid]["E_smd_eV"]) > 1e-8:
        raise ValueError("Copied SP disagrees with its qualified energy receipt")
    return sp._cds_original_ensure(gid, state_name, state, q, m, **kwargs)


def cached_pilot():
    from redox.workflows.resources import verify_gate_owner

    verify_gate_owner()
    from redox.qm.hessian_resume import resume_native
    from redox.qm.native_thermal import select_policy
    from redox.core.protocol import load_record
    from redox.core.common import apply_thermal_policy, free_energy
    from redox.qm import sp

    plan = api.read_plan(CAMPAIGN)
    api.verify_sources(plan["producer_sha256"])
    api.verify_sources(plan["input_sha256"])
    task = next(t for t in plan["tasks"] if t["key"] == KEY)
    energy, _ = api.admitted_records(CAMPAIGN, "energy")
    source = Path(energy[KEY]["state_dir"])
    state = Path(plan["raw_root"]) / "thermal" / KEY / "rest"
    state.mkdir(parents=True, exist_ok=False)
    for name in ("opt.xyz", "result.json"):
        shutil.copyfile(source / name, state / name)
    sp._cds_original_ensure = sp.ensure
    seed_electronic_record(KEY, "rest", state, 0, 1, backend="gpu")
    select_policy(state)
    report_path = CAMPAIGN / "thermal/tasks" / KEY / "status.json"
    report_path.parent.mkdir(exist_ok=False)
    owner = dict(
        host=socket.gethostname(),
        pid=os.getpid(),
        process_start=Path("/proc/self/stat").read_text().rsplit(")", 1)[1].split()[19],
        gpu=os.environ["CUDA_VISIBLE_DEVICES"],
        cpus=sorted(os.sched_getaffinity(0)),
    )
    report = dict(
        task=task,
        status="computing",
        stage="native_thermal",
        owner=owner,
        started=time.time(),
    )
    atomic_write_json(report_path, report)
    try:
        spec = json.loads((CAMPAIGN / "cached_hessian.json").read_text())
        row = dict(id=KEY, state="rest", charge=0, mult_hint=1, smiles=task["smiles"])
        resume_native(KEY, "rest", state, spec, row)
        # Preserve the original analytic-curvature RAM upper bound for dispatch.
        original = (
            Path(spec["source_hessian"]).parents[2] / "solvent_memory_execution.json"
        )
        proof = json.loads(original.read_text())
        proof.update(
            source_parent=str(original),
            scope=(
                "Original Hessian execution resource evidence; "
                "corrected CDS lineage is separate"
            ),
        )
        atomic_write_json(state / "solvent_memory_execution.json", proof)
        active = load_record(state, 0, 1)
        raw = json.loads((state / "result.json").read_text())
        qualified = apply_thermal_policy(dict(raw, e_smd_eV=active["e_smd_eV"]), state)
        G, qc = free_energy(qualified)
        if G is None or qc != "ok":
            raise ValueError("Recovered native free energy failed admission")
        api.verify_sources(plan["producer_sha256"])
        api.verify_sources(plan["input_sha256"])
        report.update(
            status="state_records_verified",
            state_dir=str(state),
            G_eV=G,
            E_smd_eV=active["e_smd_eV"],
            thermal_QC=qc,
            thermal_definition=qualified["thermal_definition"],
            finished=time.time(),
            output_sha256={
                str(p): api.digest(p)
                for p in state.rglob("*")
                if p.is_file() and p.suffix in (".json", ".xyz", ".npy")
            },
        )
        atomic_write_json(report_path, report)
        receipt = CAMPAIGN / "thermal/done" / (KEY + ".json")
        atomic_write_json(
            receipt,
            dict(
                key=KEY,
                status_path=str(report_path),
                status_sha256=api.digest(report_path),
                owner=owner,
            ),
        )
        verify_state_checkpoint(receipt)
        atomic_write_json(
            CAMPAIGN / "thermal/pilot/status.json",
            dict(
                status="thermal_pilot_passed",
                task=task,
                reused_result_receipt=str(receipt),
                receipt_sha256=api.digest(receipt),
                scope=(
                    "Same-geometry saved-Hessian recovery; all native and "
                    "worst-column gates passed"
                ),
            ),
        )
    except BaseException as error:
        report.update(
            status="failed_preserved", error=repr(error), finished=time.time()
        )
        atomic_write_json(report_path, report)
        atomic_write_json(CAMPAIGN / "thermal/pilot/status.json", report)
        raise
    finally:
        api.collect(CAMPAIGN)


def runner():
    from redox.qm import projected_thermal as projected, sp
    from redox.qm.hessian_resume import annotate_finer_hessian

    module = api.load_module(
        "CDF_finer_memory_parent",
        ROOT / "scripts/pipeline/resume_bounded_polymer_thermal.py",
    )
    module.CAMPAIGN = CAMPAIGN
    previous = projected.ensure_projected_thermal

    def seal(path, geometry):
        annotate_finer_hessian(path, geometry)
        return previous(path, geometry)

    projected.ensure_projected_thermal = seal
    try:
        backend = module.configure().runner()
    finally:
        projected.ensure_projected_thermal = previous
    backend.SCRIPT = SCRIPT
    if not hasattr(sp, "_cds_original_ensure"):
        sp._cds_original_ensure = sp.ensure
    sp.ensure = seed_electronic_record
    return backend


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "action",
        choices=(
            "prepare",
            "cached-pilot",
            "worker",
            "launch",
            "coordinate",
            "collect",
        ),
    )
    parser.add_argument("--campaign", type=Path, default=CAMPAIGN)
    parser.add_argument("--stage", choices=("thermal",), default="thermal")
    parser.add_argument("--pilot", action="store_true")
    parser.add_argument("--max-hours", type=float, default=48)
    args = parser.parse_args()
    if args.campaign.resolve() != CAMPAIGN:
        parser.error("Use the frozen CDS continuation campaign")
    if args.action == "prepare":
        prepare()
    elif args.action == "cached-pilot":
        cached_pilot()
    elif args.action == "collect":
        api.collect(CAMPAIGN)
    elif args.action == "coordinate":
        pilot = json.loads((CAMPAIGN / "thermal/pilot/status.json").read_text())
        if pilot["status"] != "thermal_pilot_passed":
            raise RuntimeError("Cached native pilot must pass before queue expansion")
        runner().coordinate(CAMPAIGN, args.max_hours)
    elif args.action == "worker":
        from redox.qm.cds_finite import finer_cds_hessian

        backend = runner()
        with finer_cds_hessian(1e-5):
            backend.worker(CAMPAIGN, "thermal", args.pilot)
    else:
        runner().launch(CAMPAIGN, "thermal", args.pilot)


if __name__ == "__main__":
    main()
