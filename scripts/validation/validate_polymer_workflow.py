#!/usr/bin/env python
"""Independently reconstruct the finite polymer screen from qualified raw states.

This does not use the production nonadditivity/conditional-grafting arithmetic.
Missing thermal states stay missing. Completion requires the selected geometry's
own native Hessian, implementation replay evidence and field-specific QC.
"""

import argparse
from collections import Counter, defaultdict
import hashlib
import json
import math
from pathlib import Path
import time

from rdkit import Chem

from redox.core.common import ROOT, EV_KJ
from redox.core.io import atomic_write_json
from redox.core.protocol import load_record
from redox.qm.native_thermal import load_native_thermal, cache_path
from redox.workflows import polymer_campaign as api
from redox.workflows.checkpoints import verify_state_checkpoint

_CPU_BINDING = None
VALIDATOR_SOURCE_BYTES = Path(__file__).read_bytes()
VALIDATOR_SOURCE_SHA256 = hashlib.sha256(VALIDATOR_SOURCE_BYTES).hexdigest()


def bind_validation_cpu():
    """Retain this validator's allocation across calls and its owned CPU client.

    A cached thermochemistry subprocess inherits the parent's assigned CPU. It
    must not be mistaken for a newly conflicting worker on a subsequent review.
    """
    import os
    from redox.workflows.resources import (
        allocation_guard,
        host_capacity,
        cpu_allocation,
    )

    global _CPU_BINDING
    with allocation_guard():
        host_capacity(1)
        if _CPU_BINDING is not None and _CPU_BINDING[0] == os.getpid():
            if tuple(sorted(os.sched_getaffinity(0))) != _CPU_BINDING[1]:
                raise RuntimeError("Validator CPU affinity changed after allocation")
            return
        cpus = cpu_allocation(1)
        os.sched_setaffinity(0, cpus)
        _CPU_BINDING = (os.getpid(), tuple(sorted(cpus)))


def finite(value):
    return type(value) in (int, float) and math.isfinite(value)


def molecular_counts(smiles):
    molecule = Chem.AddHs(Chem.MolFromSmiles(smiles))
    return Counter(
        atom.GetSymbol() for atom in molecule.GetAtoms()
    ), Chem.GetFormalCharge(molecule)


def balance(terms):
    atoms = defaultdict(int)
    charge = 0
    for coefficient, smiles in terms:
        elements, q = molecular_counts(smiles)
        for symbol, count in elements.items():
            atoms[symbol] += coefficient * count
        charge += coefficient * q
    if any(atoms.values()) or charge:
        raise ValueError(f"Unbalanced reaction/cycle: {dict(atoms)}, charge {charge}")


def references(campaign, plan):
    terms = json.loads((campaign / "reaction_terms.json").read_text())
    api.verify_sources(terms["source_sha256"])
    components = {}
    for name in terms["source_sha256"]:
        path = Path(name)
        if path.parent.name != "done":
            continue
        report = verify_state_checkpoint(path)
        key = report["task"]["key"]
        smiles = report["task"].get("smiles") or report["task"].get(
            "manifest_row", {}
        ).get("smiles")
        components[key] = dict(G_eV=report["G_eV"], smiles=smiles)
    base = terms["base"]
    precursor = [
        row
        for key, row in components.items()
        if Chem.MolToSmiles(Chem.MolFromSmiles(row["smiles"]))
        == terms["precursor_smiles"]
    ]
    if len(precursor) != 1:
        raise ValueError("Unique exact precursor reference missing")
    T = terms["temperature_K"]
    if (
        T != plan["thermochemistry"]["temperature_K"]
        or T != plan["grafting"]["temperature_K"]
    ):
        raise ValueError("Temperature mismatch")
    reagents = dict(
        precursor=precursor[0],
        base=components[base + "_neutral"],
        protonated=components[base + "_protonated"],
        pair=components[base + "_chloride_pair"],
        chloride=components["chloride"],
    )
    free = (
        reagents["protonated"]["G_eV"]
        + reagents["chloride"]["G_eV"]
        - reagents["base"]["G_eV"]
        - reagents["precursor"]["G_eV"]
    )
    standard = -0.00831446261815324 * T * math.log(0.08205736608096 * T)
    paired = (
        reagents["pair"]["G_eV"]
        - reagents["base"]["G_eV"]
        - reagents["precursor"]["G_eV"]
        + standard / EV_KJ
    )
    reconstructed = {"free_ions": free, "contact_pair": paired}
    for name, value in reconstructed.items():
        if abs(value - terms["constant_eV"][name]) > 1e-8:
            raise ValueError("Shared reagent constant does not reconstruct")
    return reconstructed, reagents


def thermal_check(report):
    state = Path(report["state_dir"])
    energy = load_record(state, 0, 1)
    thermal = load_native_thermal(state, 0, 1)
    if energy is None or thermal is None:
        raise ValueError(
            "Backing active energy or native thermochemistry fails admission"
        )
    G = energy["e_smd_eV"] + thermal["g_thermal_total_eV"]
    if abs(energy["e_smd_eV"] - report["E_smd_eV"]) > 1e-8:
        raise ValueError("Reported thermal-state electronic energy differs from its SP")
    if not finite(G) or abs(G - report["G_eV"]) > 1e-7:
        raise ValueError("Reported G is inconsistent with raw E and thermal correction")
    path, _ = cache_path(state, 0, 1)
    raw = json.loads(path.read_text())
    execution = raw.get("execution_bounded_pcm")
    if (
        not execution
        or execution.get("scientific_model_changed") is not False
        or execution.get("surface_grid_changed") is not False
    ):
        raise ValueError("Bounded PCM execution evidence missing or wrong scope")
    module = ROOT / "src/redox/qm/pcm_bounded.py"
    driver = ROOT / "scripts/pipeline/resume_bounded_polymer_thermal.py"
    if execution.get("driver_sha256") != api.digest(driver):
        raise ValueError("Bounded derivative driver or execution mode differs")
    if (
        execution.get("mode")
        == "cached_fused_PCM_Hessian_with_qualified_CDS_replacement"
    ):
        from redox.qm.hessian_resume import verify_reconstruction_record

        verify_reconstruction_record(raw, state / "opt.xyz")
    elif execution.get("mode") != "fused_diagonal_and_CPU_primitive_surface_blocks":
        raise ValueError("Unreviewed bounded derivative execution mode")
    derivative = raw.get("execution_derivatives") or {}
    for field, target in [
        ("execution_helper_sha256", ROOT / "src/redox/qm/hessian_resume.py"),
        ("finite_CDS_module_sha256", ROOT / "src/redox/qm/cds_finite.py"),
    ]:
        if field in derivative and derivative[field] != api.digest(target):
            raise ValueError("Finer-CDS execution source changed")
    if execution["adapter_sha256"] != api.digest(module) or execution[
        "source_geometry_sha256"
    ] != api.digest(state / "opt.xyz"):
        raise ValueError("Bounded execution source or geometry mismatch")
    api.verify_sources(execution["replay_evidence_sha256"])
    multiplicities = set()
    for name in execution["replay_evidence_sha256"]:
        replay = json.loads(Path(name).read_text())
        multiplicities.add(replay.get("mult"))
        if (
            replay["status"] != "full_native_replay_passed"
            or not replay["native_gates_passed"]
        ):
            raise ValueError("Implementation replay failed")
        api.verify_sources(replay["evidence_sha256"])
        if any(
            not finite(value) or not 0 <= value <= 1e-6
            for value in replay["errors"].values()
        ):
            raise ValueError("Implementation replay exceeds its fixed tolerance")
    if multiplicities != {1, 2}:
        raise ValueError("Both singlet and radical full-Hessian replays are required")
    api.verify_sources(
        {
            str(ROOT / name): value
            for name, value in execution["target_evidence_sha256"].items()
        }
    )
    corrections = {
        cutoff: energy["e_smd_eV"] + correction + thermal["g_spin_partition_eV"]
        for cutoff, correction in thermal["nuclear_by_cutoff_eV"].items()
    }
    return dict(
        G_eV=G,
        cutoff_G_eV=corrections,
        native_record=str(path),
        native_record_sha256=api.digest(path),
        minimum_frequency_cm=min(thermal["frequencies_cm"]),
        gradient_max_Ha_Bohr=thermal["gradient_max"],
        maximum_asymmetry=thermal["maximum_asymmetry"],
        derivative_checks=thermal["derivative_checks"],
    )


def validate(campaign, out, require_complete=False):
    bind_validation_cpu()
    plan = api.read_plan(campaign)
    api.verify_sources(plan["producer_sha256"])
    api.verify_sources(plan["input_sha256"])
    energy, invalid = api.admitted_records(campaign, "energy")
    thermal, bad = api.admitted_records(campaign, "thermal")
    invalid += bad
    extra_energy = (
        ROOT
        / "results/polymer_models/campaigns/2026-10-21_spacing_singlepoint_repair"
        / "energy/done/occ_none__pose1.json"
    )
    if extra_energy.exists():
        record = verify_state_checkpoint(extra_energy)
        if record.get("physical_QC") != "passed":
            raise ValueError("Repaired energy lacks physical QC")
        expected = next(t for t in plan["tasks"] if t["key"] == record["task"]["key"])
        if any(record["task"].get(key) != value for key, value in expected.items()):
            raise ValueError(
                "Repaired energy differs from the frozen molecular identity"
            )
        energy[record["task"]["key"]] = record
    extra_thermal = (
        ROOT
        / "results/polymer_models/campaigns/2026-10-21_spacing_extra_thermal"
        / "thermal/done"
    )
    campaigns = [extra_thermal.parent.parent]
    registry = ROOT / "results/polymer_models/thermal_continuations.json"
    if registry.exists():
        campaigns += [
            ROOT / path for path in json.loads(registry.read_text())["campaigns"]
        ]
    for extra_campaign in campaigns:
        if not (extra_campaign / "plan.json").exists():
            continue
        extra_plan = api.read_plan(extra_campaign)
        api.verify_sources(extra_plan["producer_sha256"])
        api.verify_sources(extra_plan["input_sha256"])
        extra_records, extra_invalid = api.admitted_records(extra_campaign, "thermal")
        invalid += extra_invalid
        expected_tasks = {t["key"]: t for t in plan["tasks"]}
        for key, report in extra_records.items():
            expected = expected_tasks.get(key)
            if expected is None or any(
                report["task"].get(k) != v for k, v in expected.items()
            ):
                raise ValueError(
                    "Extra thermal receipt differs from the original exact task"
                )
            thermal[key] = report
    if invalid:
        raise ValueError(f"Invalid receipts: {invalid}")
    constants, reagents = references(campaign, plan)
    audit = api.load_module(
        "independent_polymer_geometry",
        ROOT / "scripts/validation/audit_polymer_geometries.py",
    )
    selected = {}
    for record in energy.values():
        backing = load_record(Path(record["state_dir"]), 0, 1)
        if backing is None or abs(backing["e_smd_eV"] - record["E_smd_eV"]) > 1e-8:
            raise ValueError("Electronic report differs from its admitted raw SP")
        graph = record["task"]["model_key"]
        if graph not in selected or record["E_smd_eV"] < selected[graph]["E_smd_eV"]:
            selected[graph] = record
    cases = []
    for case in plan["graph_cases"]:
        record = selected.get(case["key"])
        row = dict(
            case,
            E_eV=None,
            G_eV=None,
            cutoff_G_eV=None,
            selected_conformer=None,
            geometry_QC="pending",
            native_QC="pending",
        )
        if record:
            review = audit.audit(record["task"], Path(record["state_dir"]) / "opt.xyz")
            if review["issues"]:
                raise ValueError(
                    f"Saved geometry failed independent review: {review['issues']}"
                )
            row.update(
                E_eV=record["E_smd_eV"],
                selected_conformer=record["task"]["key"],
                geometry_QC="passed",
            )
            free = thermal.get(record["task"]["key"])
            if free:
                if abs(free["E_smd_eV"] - record["E_smd_eV"]) > 1e-8 or api.digest(
                    Path(free["state_dir"]) / "opt.xyz"
                ) != api.digest(Path(record["state_dir"]) / "opt.xyz"):
                    raise ValueError(
                        "Selected electronic and thermal geometries differ"
                    )
                checked = thermal_check(free)
                row.update(checked, native_QC="passed")
        cases.append(row)
    lookup = {case["key"]: case for case in cases}
    cycles = []
    for definition in plan["cycles"]:
        terms = [
            (coefficient, lookup[key]["smiles"])
            for key, coefficient in definition["coefficients"].items()
        ]
        balance(terms)
        i, j = definition["sites"]
        names = ("occ_none", f"occ_{i}", f"occ_{j}", f"occ_{i}_{j}")
        values = [lookup[name] for name in names]
        E = [row["E_eV"] for row in values]
        G = [row["G_eV"] for row in values]
        delta_E = (
            (E[3] + E[0] - E[1] - E[2]) * EV_KJ if all(finite(v) for v in E) else None
        )
        delta_G = (
            (G[3] + G[0] - G[1] - G[2]) * EV_KJ if all(finite(v) for v in G) else None
        )
        route_values = {}
        for name, constant in constants.items():
            if delta_G is None:
                route_values[name] = None
                continue
            first_j = (G[2] - G[0] + constant) * EV_KJ
            second_j = (G[3] - G[1] + constant) * EV_KJ
            if abs(second_j - first_j - delta_G) > 1e-6:
                raise ValueError(
                    "Conditional grafting cycle failed independent closure"
                )
            route_values[name] = dict(
                first_i=(G[1] - G[0] + constant) * EV_KJ,
                first_j=first_j,
                second_j_given_i=second_j,
                second_i_given_j=(G[3] - G[2] + constant) * EV_KJ,
            )
        for start, finish in ((names[0], names[1]), (names[1], names[3])):
            balance(
                [
                    (1, lookup[start]["smiles"]),
                    (1, reagents["precursor"]["smiles"]),
                    (1, reagents["base"]["smiles"]),
                    (-1, lookup[finish]["smiles"]),
                    (-1, reagents["protonated"]["smiles"]),
                    (-1, reagents["chloride"]["smiles"]),
                ]
            )
        cutoff = None
        if delta_G is not None:
            cutoff = {
                tag: (
                    values[3]["cutoff_G_eV"][tag]
                    + values[0]["cutoff_G_eV"][tag]
                    - values[1]["cutoff_G_eV"][tag]
                    - values[2]["cutoff_G_eV"][tag]
                )
                * EV_KJ
                for tag in ("25", "50", "100")
            }
        cycles.append(
            dict(
                definition,
                electronic_deltaDeltaE_kJ_mol=delta_E,
                free_energy_deltaDeltaG_kJ_mol=delta_G,
                conditional_grafting_kJ_mol=route_values,
                rotor_cutoff_sensitivity_kJ_mol=cutoff,
                atoms_and_charge_balanced=True,
            )
        )
    original = json.loads(
        (ROOT / "results/polymer_models/interaction_status.json").read_text()
    )
    ready_G = sum(row["G_eV"] is not None for row in cases)
    complete = (
        original["stationary_electronic_cases"] == original["expected_DFT_cases"]
        and len(energy) == len(plan["tasks"])
        and ready_G == len(cases)
    )
    source_directory = Path(out).parent / "validator_sources"
    source_directory.mkdir(exist_ok=True)
    source_snapshot = source_directory / (VALIDATOR_SOURCE_SHA256 + ".py")
    if not source_snapshot.exists():
        source_snapshot.write_bytes(VALIDATOR_SOURCE_BYTES)
    if api.digest(source_snapshot) != VALIDATOR_SOURCE_SHA256:
        raise ValueError("Validator source snapshot changed")
    result = dict(
        status="finite_workflow_validated"
        if complete
        else "incomplete_with_explicit_pending_fields",
        checked_at=time.time(),
        campaign=str(campaign),
        original_checked=original["stationary_electronic_cases"],
        original_expected=original["expected_DFT_cases"],
        qualified_spacing_conformers=len(energy),
        expected_spacing_conformers=len(plan["tasks"]),
        qualified_native_G=ready_G,
        expected_graphs=len(cases),
        cases=cases,
        cycles=cycles,
        source_plan_sha256=api.digest(campaign / "plan.json"),
        selected_geometry_G_pairing_required=True,
        no_missing_thermal_substitution=True,
        limits=[
            "Finite neutral method-control screen; no experimental synthesis admission",
            "Single-conformer thermochemistry; no polymer occupancy probabilities",
            "Nonadditivity includes backbone/solvent response; no pure-core binding",
            "Rotor-cutoff variation is sensitivity, not a confidence interval",
        ],
        validator_sha256=VALIDATOR_SOURCE_SHA256,
        validator_source_snapshot=str(source_snapshot),
    )
    atomic_write_json(out, result)
    print(
        json.dumps(
            {
                key: result[key]
                for key in (
                    "status",
                    "original_checked",
                    "qualified_spacing_conformers",
                    "qualified_native_G",
                )
            }
        )
    )
    if require_complete and not complete:
        raise RuntimeError("Workflow is not complete; pending fields remain explicit")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--campaign",
        type=Path,
        default=ROOT
        / "results/polymer_models/campaigns/2026-10-21_spacing_bounded_pcm",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=ROOT / "results/polymer_models/workflow_validation.json",
    )
    parser.add_argument("--require-complete", action="store_true")
    args = parser.parse_args()
    validate(args.campaign.resolve(), args.out.resolve(), args.require_complete)
