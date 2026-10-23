#!/usr/bin/env python
"""Freeze final CREST conformers as reproducible DFT re-ranking inputs.

Running searches, MD frames and intermediate ensemble files never enter this
queue. Every selected geometry is bound to its exact occupancy/stereo graph.
The initial ordered batches contain 10 then 20 starts, retaining both previously
checked DFT poses where available before adding final CREST energy-ranked starts.
"""
import argparse
import json
import math
import os
from pathlib import Path
import shutil
import time

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api

SAMPLING = ROOT / "results/validation/conformer_sampling_2026-10-23"
FRAGMENT = "[*:1]N1C(=O)c2ccccc2C1=O"


def plain_xyz_energy(frame):
    """CREST plain-XYZ comments carry Hartree, unlike extended XYZ defaults."""
    comment = frame.splitlines()[1].strip()
    try:
        value = float(comment)
    except ValueError as error:
        raise ValueError("Final plain-XYZ comment must contain only an energy in Hartree") from error
    if not math.isfinite(value):
        raise ValueError("Nonfinite CREST ranking energy")
    return value


def configuration(case):
    tags = case["backbone_stereo_tags"]
    if tags == ["CW"] * 4:
        return "all_same"
    if tags == ["CW", "CCW", "CW", "CCW"]:
        return "alternating"
    raise ValueError("Unregistered configuration in the current six-graph panel")


def prepare_request(root):
    root.mkdir(parents=True, exist_ok=False)
    initial = SAMPLING / "inputs/plan.json"
    wave = SAMPLING / "remaining_wave_plan.json"
    tasks = json.loads(initial.read_text())["tasks"]
    tasks += [row for row in json.loads(wave.read_text())["tasks"] if "parallel_repeat" in row["key"]]
    if len(tasks) != 12:
        raise ValueError("Expected twelve primary graph/configuration searches")
    sources = {str(path): api.digest(path) for path in (initial, wave)}
    baseline = ROOT / "results/polymer_models/campaigns/2026-10-21_spacing_cds_recovery"
    reference_energy, invalid = api.admitted_records(baseline, "energy")
    if invalid or len(reference_energy) != 12:
        raise ValueError("Previously checked reference poses must remain available")
    cases = []
    for task in tasks:
        path = Path(task["case"])
        case = json.loads(path.read_text())
        label = configuration(case)
        occupied = "_".join(map(str, case["occupied_sites"])) or "none"
        key = label + "__occ_" + occupied
        sources[str(path)] = api.digest(path)
        references = []
        if label == "all_same":
            for ref_key, record in sorted(reference_energy.items()):
                if record["task"]["smiles"] != case["smiles"]:
                    continue
                receipt = baseline / "energy/done" / (ref_key + ".json")
                xyz = Path(record["state_dir"]) / "opt.xyz"
                references.append(dict(key=ref_key, path=str(xyz), sha256=api.digest(xyz), receipt=str(receipt)))
                sources.update({str(p): api.digest(p) for p in (receipt, xyz)})
            if len(references) != 2:
                raise ValueError("Both earlier exact-graph DFT poses must be retained")
        cases.append(dict(key=key, configuration=label, case_path=str(path), source_search=task["raw_out"],
            mandatory_references=references,
            fragment=FRAGMENT, motif_id="exp_phthalimide_control", role="interaction_method_control; not synthesis admission"))
    if len({row["key"] for row in cases}) != 12:
        raise ValueError("Duplicate graph/configuration request")
    request = dict(schema=1, cases=cases, ordered_DFT_batches=[10, 20], maximum_initial_starts=20,
        temperature_K=298.15, solvent="acetonitrile", source_sha256=sources,
        source_roles="Twelve primary searches; earlier slow reference searches remain separate sensitivity inputs",
        raw_root=str(ROOT / "calcs/polymer_models" / root.name), created=time.time(),
        experimental_loading_known=False, synthesis_admission=False,
        scope="Neutral phthalimide occupancy/tacticity method control; no experimental population weights")
    atomic_write_json(root / "request.json", request)
    return request


def source_records(entry):
    """Verify final source completion and every frame before selecting starts."""
    from redox.build.geometry import geometry_stereo_status
    from redox.validation.protonation import hydrogen_topology_status

    raw = Path(entry["source_search"])
    status_path = raw / "status.json"
    status = json.loads(status_path.read_text())
    if status["status"] != "low_level_ensemble_complete_DFT_pending":
        return None
    if not status.get("finished"):
        raise ValueError("A final search lacks its terminal source record")
    producer = str(ROOT / "scripts/pipeline/sample_polymer_conformers.py")
    if producer not in status.get("source_sha256", {}) or str(Path(entry["case_path"]).resolve()) not in status["source_sha256"]:
        raise ValueError("Final source lacks its sampler/case provenance")
    api.verify_sources(status["source_sha256"])
    argv = status.get("argv", [])
    for option, value in (("--alpb", "acetonitrile"), ("--chrg", "0"), ("--uhf", "0")):
        if option not in argv or argv.index(option) + 1 >= len(argv) or argv[argv.index(option) + 1] != value:
            raise ValueError("Final source medium/charge/spin differs")
    if "--gfn2" not in argv:
        raise ValueError("Final source uses a different low-level method")
    case = json.loads(Path(entry["case_path"]).read_text())
    if status["case"] != case or case.get("charge") != 0 or case.get("mult") != 1:
        raise ValueError("Source identity differs or needs an unqualified ionic/spin model")
    from redox.build.polymer_cases import build_case
    rebuilt = build_case(entry["fragment"], case["repeat_units"], case["reactive_sites"], case["occupied_sites"], case["backbone_stereo_tags"])
    if rebuilt["smiles"] != case["smiles"] or rebuilt["backbone_stereo_fingerprint"] != case["backbone_stereo_fingerprint"]:
        raise ValueError("Source is not the requested open-chain occupancy/stereo graph")
    metadata_path = raw / "conformers.json"
    metadata = json.loads(metadata_path.read_text())
    ensemble = raw / "crest_conformers.xyz"
    if api.digest(ensemble) != metadata["ensemble_sha256"]:
        raise ValueError("Final CREST ensemble changed")
    reader = api.load_module("refinement_xyz_frames", ROOT / "scripts/pipeline/sample_polymer_crest.py")
    frames = list(reader.xyz_frames(ensemble))
    if len(frames) != len(metadata["records"]) or len(frames) != status["conformers"]:
        raise ValueError("Final ensemble/frame inventory differs")
    valid, rejected = [], []
    import hashlib
    for index, (frame, row) in enumerate(zip(frames, metadata["records"])):
        path = Path(row["path"])
        if row["index"] != index or not path.resolve().is_relative_to(raw.resolve()):
            raise ValueError("Conformer index/path binding differs")
        if hashlib.sha256(frame.encode()).hexdigest() != row["sha256"] or api.digest(path) != row["sha256"]:
            raise ValueError("Conformer geometry differs from its final source frame")
        graph = hydrogen_topology_status(case["smiles"], path)
        stereo = geometry_stereo_status(case["smiles"], path)
        good = graph == "ok" and stereo in ("ok", "not_applicable")
        if good != row["eligible_for_DFT_reranking"]:
            raise ValueError("Independent graph/stereo replay differs")
        record = dict(index=index, path=str(path), sha256=row["sha256"], GFN2_energy_Ha=plain_xyz_energy(frame), graph=graph, stereo=stereo)
        (valid if good else rejected).append(record)
    if not valid:
        raise ValueError("Final search contains no eligible exact-graph conformer")
    sources = dict(status["source_sha256"])
    sources.update({str(path): api.digest(path) for path in (status_path, metadata_path, ensemble, Path(entry["case_path"]))})
    return case, sorted(valid, key=lambda row: (row["GFN2_energy_Ha"], row["index"])), rejected, sources


def stage_case(root, request, entry):
    campaign = root / entry["key"]
    if (campaign / "plan.json").exists():
        plan = api.read_plan(campaign)
        api.verify_sources(plan["producer_sha256"])
        api.verify_sources(plan["input_sha256"])
        return dict(key=entry["key"], status="DFT_inputs_frozen", starts=len(plan["tasks"]))
    source = source_records(entry)
    if source is None:
        state = json.loads((Path(entry["source_search"]) / "status.json").read_text())["status"]
        return dict(key=entry["key"], status="source_search_failed_preserved" if state.startswith("failed") else "waiting_for_final_CREST", source_status=state, starts=0)
    case, ranked, rejected, input_hashes = source
    if campaign.exists():
        raise RuntimeError("Unfinished preparation must be reviewed; preserve its files")
    campaign.mkdir()
    from rdkit import Chem
    from rdkit.Chem import rdDetermineBonds
    from ase.io import read
    from pyscf import gto
    from redox.core.protocol import ACTIVE_SP
    from redox.build.geometry import _skeleton
    from redox.build.polymer_contacts import pi_core_groups
    audit = api.load_module("refinement_input_geometry", ROOT / "scripts/validation/audit_polymer_geometries.py")
    references = []
    from redox.workflows.checkpoints import verify_state_checkpoint
    for reference in entry.get("mandatory_references", []):
        record = verify_state_checkpoint(reference["receipt"])
        if record["task"]["smiles"] != case["smiles"] or api.digest(reference["path"]) != reference["sha256"]:
            raise ValueError("Previously checked reference identity/geometry differs")
        references.append(dict(reference, index=None, GFN2_energy_Ha=None, source_role="previous_checked_DFT_pose"))
        input_hashes[reference["receipt"]] = api.digest(reference["receipt"])
    selected = references + [dict(row, source_role="final_CREST_conformer") for row in ranked[:request["maximum_initial_starts"] - len(references)]]
    declared = Chem.AddHs(Chem.MolFromSmiles(case["smiles"]))
    if declared.GetNumAtoms() != case["atoms"]:
        raise ValueError("Case atom count differs")
    template_groups = pi_core_groups(declared, entry["fragment"])
    if len(template_groups) != len(case["occupied_sites"]):
        raise ValueError("Conjugated motif count differs from occupancy")
    atom_labels = None
    tasks = []
    for rank, row in enumerate(selected, 1):
        key = entry["key"] + f"__conf{rank:03d}"
        folder = campaign / "live_inputs" / key
        folder.mkdir(parents=True)
        xyz = folder / "input.xyz"
        shutil.copyfile(row["path"], xyz)
        atoms = read(xyz)
        labels = atoms.get_chemical_symbols()
        if atom_labels is not None and labels != atom_labels:
            raise ValueError("CREST changed the atom ordering between conformers")
        atom_labels = labels
        observed = Chem.MolFromXYZFile(str(xyz))
        rdDetermineBonds.DetermineConnectivity(observed, charge=0)
        mapping = _skeleton(observed).GetSubstructMatch(_skeleton(declared))
        if len(mapping) != declared.GetNumAtoms():
            raise ValueError("Cannot map every declared atom to the saved XYZ")
        groups = [[mapping[i] for i in group] for group in template_groups]
        basis_count = gto.M(atom=list(zip(labels, atoms.positions)), basis=ACTIVE_SP["basis"], charge=0, spin=0, unit="Angstrom", verbose=0).nao_nr()
        task = dict(case, key=key, model_key=entry["key"], candidate_id=entry["motif_id"],
            model_type="covalent_oligomer", repeats=case["repeat_units"], stereoisomer=entry["configuration"],
            target="previous_checked_DFT_reference" if row["source_role"] == "previous_checked_DFT_pose" else f"final_CREST_index_{row['index']}", fragment=entry["fragment"], pi_groups=groups,
            active_basis_functions=basis_count, mult=1, role=entry["role"], search_rank=rank,
            source_conformer_index=row["index"], source_GFN2_energy_Ha=row["GFN2_energy_Ha"])
        task["starting_geometry_source_role"] = row["source_role"]
        task["source_reference_key"] = row.get("key")
        review = audit.audit(task, xyz)
        if review["issues"]:
            raise ValueError("Selected conformer failed physical input QC: " + str(review["issues"]))
        atomic_write_json(folder / "geometry_review.json", review)
        hashes = {str(path): api.digest(path) for path in (xyz, Path(row["path"]), folder / "geometry_review.json")}
        atomic_write_json(folder / "task.json", dict(task, input_xyz=str(xyz), input_sha256=hashes))
        input_hashes.update(hashes)
        input_hashes[str(folder / "task.json")] = api.digest(folder / "task.json")
        tasks.append(task)
    for stage in ("energy", "thermal"):
        for name in ("tasks", "done", "claims", "failures", "workers", "pilot", "launches"):
            (campaign / stage / name).mkdir(parents=True)
    # The checked memory-replay controls remain source-bound inputs. They are
    # implementation controls; each new graph still requires resource/physical QC.
    parent = ROOT / "results/polymer_models/campaigns/2026-10-21_spacing_cds_recovery"
    shutil.copytree(parent / "controls", campaign / "controls")
    for path in (campaign / "controls").glob("*/status.json"):
        input_hashes[str(path)] = api.digest(path)
    input_hashes[str(root / "request.json")] = api.digest(root / "request.json")
    sources = list((ROOT / "src/redox").rglob("*.py")) + [Path(__file__).resolve()]
    sources += [ROOT / "scripts/pipeline" / name for name in (
        "refine_polymer_conformers.py", "run_polymer_spacing.py", "run_pi_stacking.py",
        "polymer_conformer_thermals.py", "resume_polymer_forces_adaptive.py", "resume_polymer_forces.py",
        "continue_polymer_models.py", "resume_bounded_polymer_thermal.py", "resume_polymer_thermal.py")]
    sources += [ROOT / "scripts/validation" / name for name in ("audit_polymer_geometries.py", "diagnose_saved_hessian.py")]
    plan = dict(schema=1, tasks=tasks, graph_cases=[dict(case, key=entry["key"])], cycles=[],
        source_search=entry["source_search"], source_eligible_conformers=len(ranked), rejected_conformers=rejected,
        ordered_DFT_batches=request["ordered_DFT_batches"], requested_initial_starts=20,
        raw_root=str(Path(request["raw_root"]) / entry["key"]), solvent="acetonitrile",
        quantum_settings=dict(max_opt_steps=400, geometry_SCF_conv_tol=1e-11),
        thermochemistry=dict(temperature_K=298.15, memory_cap_gib=24),
        minimum_large_GPU_MiB=80 * 1024, maximum_memory_pilot_MiB=80 * 1024,
        pilot_key=max(tasks, key=lambda task: task["active_basis_functions"])["key"],
        producer_sha256={str(path): api.digest(path) for path in sources}, input_sha256=input_hashes,
        created=time.time(), minimum_change_tolerance_kJ_mol=2, ensemble_change_tolerance_kJ_mol=2,
        native_G_required=True, experimental_occupancy_weights=None, synthesis_admission=False,
        interpretation="Final low-level minima are DFT warm starts, not an admitted thermodynamic ensemble")
    api.verify_sources(input_hashes)
    atomic_write_json(campaign / "plan.json", plan)
    return dict(key=entry["key"], status="DFT_inputs_frozen", starts=len(tasks))


def stage(root):
    request = json.loads((root / "request.json").read_text())
    api.verify_sources(request["source_sha256"])
    rows = []
    for entry in request["cases"]:
        try:
            rows.append(stage_case(root, request, entry))
        except (OSError, ValueError, RuntimeError) as error:
            rows.append(dict(key=entry["key"], status="preparation_requires_review_outputs_preserved", starts=0, error=repr(error)))
    result = dict(cases=rows, frozen_graphs=sum(row["starts"] > 0 for row in rows),
        staged_DFT_starts=sum(row["starts"] for row in rows), quantum_jobs_submitted=0, updated=time.time())
    atomic_write_json(root / "preparation_status.json", result)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("prepare", "stage"))
    parser.add_argument("--root", type=Path, required=True)
    args = parser.parse_args()
    from redox.workflows.resources import allocation_guard, cpu_allocation, host_capacity
    with allocation_guard():
        host_capacity(1)
        os.sched_setaffinity(0, cpu_allocation(1))
    root = args.root.resolve()
    if args.action == "prepare":
        prepare_request(root)
    result = stage(root)
    print(json.dumps({key: value for key, value in result.items() if key != "cases"}))


if __name__ == "__main__":
    main()
