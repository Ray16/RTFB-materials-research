"""Scientific state and analysis for reproducible polymer occupancy campaigns.

Geometry/SP calculation uses a thin adapter to the existing, qualified oligomer
driver. Electronic nonadditivity, conditional grafting free energies and physical
feasibility are deliberately separate fields with separate admission rules.
"""

from contextlib import contextmanager
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import shutil
import time

from redox.core.common import ROOT, EV_KJ
from redox.core.io import atomic_write_json
from redox.workflows.checkpoints import verify_state_checkpoint


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def verify_sources(mapping):
    for path, expected in mapping.items():
        if digest(path) != expected:
            raise RuntimeError(f"Source or input changed: {path}")


def load_module(name, path):
    specification = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


def read_plan(campaign):
    return json.loads((Path(campaign) / "plan.json").read_text())


def reaction_terms(campaign, motif_id, base, overrides=None, *, reference_directory=None, base_directory=None):
    """Freeze shared reagent G terms, leaving missing references explicitly blank.

    The contact-pair reaction has delta-n=-1; separated ions have delta-n=0.
    Missing reagent data blocks absolute conditional grafting G, while a balanced
    neighbor-penalty cycle can still be evaluated once the polymer G values exist.
    """
    from rdkit import Chem
    from redox.build.grafting import precursor

    inventory_path = ROOT / "results/grafting/reactions.json"
    inventory = json.loads(inventory_path.read_text()) if inventory_path.exists() else {"reactions": []}
    case_spec = json.loads((Path(campaign) / "cases.json").read_text())
    exact_precursor, exact_class = precursor(case_spec["fragment"])
    exact_smiles = Chem.MolToSmiles(exact_precursor)
    matches = [row for row in inventory["reactions"] if row["id"] == motif_id]
    if matches:
        reaction = matches[0]
        precursor_id = reaction["precursor_id"]
        precursor_smiles = reaction["precursor_smiles"]
        reaction_class = reaction["reaction_class"]
        if precursor_smiles != exact_smiles or reaction_class != exact_class:
            raise ValueError("Registered route differs from the exact configured motif")
    else:
        reaction_class = exact_class
        precursor_smiles = exact_smiles
        precursor_id = None
    if reaction_class != "NH_OH_substitution":
        raise ValueError("This backend requires a neutral N-H/O-H grafting route")
    reference = Path(reference_directory) if reference_directory else ROOT / "results/grafting/references/molecular"
    bases = Path(base_directory) if base_directory else ROOT / "results/grafting/references/bases"
    base_smiles = json.loads((ROOT / "config/grafting_references.json").read_text())["base_smiles"]
    if f"{base}_neutral" not in base_smiles:
        raise ValueError("Select a qualified DBU/DIPEA base reference")
    expected_smiles = {
        "precursor": precursor_smiles,
        "chloride": "[Cl-]",
        "base": base_smiles[f"{base}_neutral"],
        "protonated": base_smiles[f"{base}_protonated"],
        "pair": base_smiles[f"{base}_chloride_pair"],
    }
    paths = {
        "precursor": reference / f"{precursor_id}.json" if precursor_id else None,
        "chloride": reference / "chloride.json",
        "base": bases / f"{base}_neutral.json",
        "protonated": bases / f"{base}_protonated.json",
        "pair": bases / f"{base}_chloride_pair.json",
    }
    paths.update({name: Path(path) for name, path in (overrides or {}).items()})
    if set(paths) != {"precursor", "chloride", "base", "protonated", "pair"}:
        raise ValueError("Unknown reagent override role")
    energies, sources, missing = {}, {}, {}
    for name, path in paths.items():
        try:
            if path is None:
                raise ValueError("No qualified reagent receipt supplied")
            report = verify_state_checkpoint(path)
            value = report["G_eV"]
            if type(value) not in (int, float) or not math.isfinite(value):
                raise ValueError("No finite reagent G")
            readiness = report.get("readiness")
            native_ready = (
                isinstance(readiness, list)
                and readiness
                and all(row["status"] == "ready" for row in readiness)
            ) or (
                report.get("thermal_QC") == "ok"
                and report.get("thermal_definition") == "native_solution_dft_qrrho_spin"
            )
            if not native_ready and readiness != "qualified_monatomic_reference":
                raise ValueError("Native or monatomic thermochemistry QC is required")
            task = report["task"]
            smiles = task.get("smiles") or task.get("manifest_row", {}).get("smiles")
            molecule = Chem.MolFromSmiles(smiles) if smiles else None
            if molecule is None or Chem.MolToSmiles(molecule) != Chem.MolToSmiles(
                Chem.MolFromSmiles(expected_smiles[name])
            ):
                raise ValueError(f"Wrong molecular identity for reagent role {name}")
            if readiness != "qualified_monatomic_reference":
                policy = json.loads(
                    (Path(report["state_dir"]) / "thermal_policy.json").read_text()
                )
                protocol = policy["native_protocol"]
                if (
                    protocol["temperature_K"] != 298.15
                    or protocol["solvent"] != "acetonitrile"
                    or protocol["solvent_model"] != "SMD"
                ):
                    raise ValueError("Reagent temperature/solvent protocol mismatch")
            energies[name] = value
            sources[str(path)] = digest(path)
            sources.update(report["output_sha256"])
        except (OSError, ValueError, KeyError) as error:
            missing[name] = repr(error)
    temperature = 298.15
    concentration = (
        -0.00831446261815324 * temperature * math.log(0.08205736608096 * temperature)
    )
    free_ready = all(
        name in energies for name in ("protonated", "chloride", "base", "precursor")
    )
    pair_ready = all(name in energies for name in ("pair", "base", "precursor"))
    terms = {
        "free_ions": energies["protonated"]
        + energies["chloride"]
        - energies["base"]
        - energies["precursor"]
        if free_ready
        else None,
        "contact_pair": energies["pair"]
        - energies["base"]
        - energies["precursor"]
        + concentration / EV_KJ
        if pair_ready
        else None,
    }
    snapshot = dict(
        status="qualified_shared_terms" if not missing else "pending_reagent_QC",
        base=base,
        temperature_K=temperature,
        route_class=reaction_class,
        precursor_smiles=precursor_smiles,
        constant_eV=terms,
        source_sha256=sources,
        missing_reagents=missing,
        standard_state=(
            "1 M molecular model with explicit contact-pair delta-n correction"
        ),
        scope=(
            "Shared reagents only; polymer G requires qualified native thermochemistry"
        ),
    )
    atomic_write_json(Path(campaign) / "reaction_terms.json", snapshot)
    return snapshot


def freeze(campaign):
    campaign = Path(campaign).resolve()
    if (campaign / "plan.json").exists():
        raise RuntimeError(
            "Preserve the existing frozen campaign; use a new run for changes"
        )
    preparation = json.loads((campaign / "preparation.json").read_text())
    if not set(preparation["grafting"]["reaction_models"]) <= {
        "free_ions",
        "contact_pair",
    }:
        raise ValueError("Unknown byproduct reaction model")
    if preparation["solvent"] != "acetonitrile":
        raise ValueError("This qualified backend requires acetonitrile")
    verify_sources(preparation["producer_sha256"])
    verify_sources(preparation["input_sha256"])
    if preparation["held_graphs"]:
        raise RuntimeError(
            "Resolve the recorded preparation holds before freezing the matched panel"
        )
    terms = reaction_terms(
        campaign,
        preparation["motif_id"],
        preparation["grafting"]["base"],
        preparation["grafting"].get("reagent_receipts"),
    )
    if preparation["thermochemistry"]["temperature_K"] != terms["temperature_K"]:
        raise ValueError(
            "Shared reagent and polymer thermochemistry temperatures must match"
        )
    sources = list((ROOT / "src/redox").rglob("*.py"))
    sources += list((ROOT / "config").glob("*.py")) + list(
        (ROOT / "config").glob("*.json")
    )
    sources += [
        ROOT / "scripts/pipeline/run_polymer_spacing.py",
        ROOT / "scripts/pipeline/prepare_polymer_spacing.py",
        ROOT / "scripts/pipeline/run_pi_stacking.py",
        ROOT / "scripts/validation/audit_polymer_geometries.py",
    ]
    plan = dict(preparation)
    plan["producer_sha256"] = {str(path): digest(path) for path in sources}
    plan["input_sha256"] = dict(preparation["input_sha256"])
    for name in ("preparation.json", "reaction_terms.json"):
        path = campaign / name
        plan["input_sha256"][str(path)] = digest(path)
    plan["frozen"] = time.time()
    if preparation["grafting"]["temperature_K"] != terms["temperature_K"]:
        raise ValueError("Grafting and reagent reference temperatures must match")
    for stage in ("energy", "thermal"):
        for directory in (
            "tasks",
            "done",
            "claims",
            "failures",
            "workers",
            "pilot",
            "launches",
        ):
            (campaign / stage / directory).mkdir(parents=True, exist_ok=True)
    atomic_write_json(campaign / "plan.json", plan)
    collect(campaign)
    return plan


def nonadditivity(empty, first, second, both):
    values = (empty, first, second, both)
    if any(
        type(value) not in (int, float) or not math.isfinite(value) for value in values
    ):
        raise ValueError("Four finite, qualified state energies are required")
    return (both + empty - first - second) * EV_KJ


def conditional_grafting(empty, first, second, both, reagent_constant):
    """Compare adding j to empty and i-occupied backbones under identical conditions."""
    if type(reagent_constant) not in (int, float) or not math.isfinite(
        reagent_constant
    ):
        raise ValueError("A finite qualified reaction constant is required")
    penalty = nonadditivity(empty, first, second, both)
    first_i = (first - empty + reagent_constant) * EV_KJ
    first_j = (second - empty + reagent_constant) * EV_KJ
    second_j = (both - first + reagent_constant) * EV_KJ
    second_i = (both - second + reagent_constant) * EV_KJ
    if abs((second_j - first_j) - penalty) > 1e-6:
        raise ValueError("Conditional grafting cycle did not close")
    return dict(
        graft_i_on_empty_kJ_mol=first_i,
        graft_j_on_empty_kJ_mol=first_j,
        graft_j_given_i_kJ_mol=second_j,
        graft_i_given_j_kJ_mol=second_i,
        neighbor_penalty_deltaDeltaG_kJ_mol=penalty,
        interpretation=(
            "Positive penalty makes grafting less favorable; it does not alone "
            "establish that the full reaction is unfavorable"
        ),
        synthesis_feasibility=(
            "Unvalidated kinetics, selectivity, accessibility and material environment"
        ),
    )


def admitted_records(campaign, stage):
    plan = read_plan(campaign)
    expected_tasks = {task["key"]: task for task in plan["tasks"]}
    records = {}
    invalid = []
    for path in (Path(campaign) / stage / "done").glob("*.json"):
        try:
            record = verify_state_checkpoint(path)
            task = record["task"]
            key = task["key"]
            expected = expected_tasks.get(key)
            if expected is None or any(
                task.get(name) != value for name, value in expected.items()
            ):
                raise ValueError("Result identity differs from the frozen task")
            field = "G_eV" if stage == "thermal" else "E_smd_eV"
            value = record.get(field)
            if type(value) not in (int, float) or not math.isfinite(value):
                raise ValueError(f"Missing finite {field}")
            if stage == "thermal" and (
                record.get("thermal_QC") != "ok"
                or record.get("thermal_definition") != "native_solution_dft_qrrho_spin"
            ):
                raise ValueError("Unqualified or mixed thermal definition")
            if stage == "energy" and record.get("physical_QC") != "passed":
                raise ValueError("Independent physical geometry QC is required")
            records[key] = record
        except (OSError, KeyError, ValueError) as error:
            invalid.append(dict(receipt=str(path), error=repr(error)))
    return records, invalid


def selected_graph_records(campaign):
    energy, invalid = admitted_records(campaign, "energy")
    selected = {}
    for record in energy.values():
        graph = record["task"]["model_key"]
        if graph not in selected or record["E_smd_eV"] < selected[graph]["E_smd_eV"]:
            selected[graph] = record
    return selected, invalid


def energy_sampling_terminal(campaign, plan):
    root = Path(campaign) / "energy"
    return all(
        (root / "done" / f"{task['key']}.json").exists()
        or (root / "failures" / f"{task['key']}.json").exists()
        for task in plan["tasks"]
    )


def collect(campaign):
    campaign = Path(campaign).resolve()
    plan = read_plan(campaign)
    selected, invalid = selected_graph_records(campaign)
    thermal, thermal_invalid = admitted_records(campaign, "thermal")
    invalid += thermal_invalid
    sampling_complete = energy_sampling_terminal(campaign, plan)
    graphs = []
    lookup = {}
    for case in plan["graph_cases"]:
        record = selected.get(case["key"])
        free = thermal.get(record["task"]["key"]) if record else None
        row = dict(
            case,
            E_eV=record["E_smd_eV"] if record and sampling_complete else None,
            G_eV=free.get("G_eV") if free and sampling_complete else None,
            selected_conformer=record["task"]["key"] if record else None,
            geometry_status="qualified_sampled_structure" if record else "pending",
            G_status="qualified_single_conformer_native_G"
            if free
            else "pending_native_thermal_QC",
        )
        graphs.append(row)
        lookup[case["key"]] = row
    terms = json.loads((campaign / "reaction_terms.json").read_text())
    verify_sources(terms["source_sha256"])
    cycles = []
    for cycle in plan["cycles"]:
        i, j = cycle["sites"]
        keys = ["occ_none", f"occ_{i}", f"occ_{j}", f"occ_{i}_{j}"]
        energy = [lookup[key]["E_eV"] for key in keys]
        free = [lookup[key]["G_eV"] for key in keys]
        e_ready = all(value is not None for value in energy)
        g_ready = all(value is not None for value in free)
        models = {
            name: conditional_grafting(*free, terms["constant_eV"][name])
            if g_ready and terms["constant_eV"][name] is not None
            else None
            for name in plan["grafting"]["reaction_models"]
        }
        cycles.append(
            dict(
                cycle,
                electronic_neighbor_penalty_deltaDeltaE_kJ_mol=nonadditivity(*energy)
                if e_ready
                else None,
                free_energy_neighbor_penalty_deltaDeltaG_kJ_mol=nonadditivity(*free)
                if g_ready
                else None,
                conditional_grafting=models,
                status="qualified_single_conformer_G_screen"
                if g_ready
                else "electronic_screen_only"
                if e_ready
                else "pending_calculations",
            )
        )
    status = dict(
        updated=time.time(),
        graph_cases=len(graphs),
        DFT_tasks=len(plan["tasks"]),
        completed_DFT=len(admitted_records(campaign, "energy")[0]),
        completed_native_G=len(thermal),
        energy_sampling_terminal=sampling_complete,
        invalid_receipts=invalid,
        active_claims={
            stage: [
                path.name
                for path in (campaign / stage / "claims").iterdir()
                if path.is_dir()
            ]
            for stage in ("energy", "thermal")
        },
        failures={
            stage: len(list((campaign / stage / "failures").glob("*.json")))
            for stage in ("energy", "thermal")
        },
        raw_root=plan["raw_root"],
        source_plan_sha256=digest(campaign / "plan.json"),
        physical_scope=(
            "Conditional molecular proxy; no experimental occupancy weights, "
            "whole-bead free energy or electron-transfer rate"
        ),
    )
    current = dict(
        status=status,
        cases=graphs,
        cycles=cycles,
        conditions=dict(
            solvent=plan["solvent"],
            grafting=plan["grafting"],
            temperature_K=terms["temperature_K"],
        ),
        conformer_selection=(
            "Lowest qualified sampled electronic-energy conformer per graph; "
            "no conformational free-energy ensemble"
        ),
        metrics_are_distinct=[
            "electronic site nonadditivity",
            "conditional grafting deltaG",
            "electronic coupling/rates remain uncomputed",
        ],
    )
    atomic_write_json(campaign / "current.json", current)
    atomic_write_json(campaign / "status.json", status)
    return current


@contextmanager
def optimizer_settings(settings, checkpoint_directory):
    """Apply run-local settings while preserving shared live producer files."""
    import redox.qm.dft as dft
    from redox.workflows.geometry_checkpoints import write_geometry

    original_dft = dft.dft_smd
    original_builder = dft._build_smd_mf
    original_optimizer = dft._geomopt_fn

    def evaluate(*args, **kwargs):
        if kwargs.get("do_opt", True):
            kwargs["max_opt_steps"] = settings["max_opt_steps"]
        return original_dft(*args, **kwargs)

    def build(*args, **kwargs):
        mean_field = original_builder(*args, **kwargs)
        if (args[1] if len(args) > 1 else kwargs.get("xc")) == "r2scan":
            mean_field.conv_tol = settings["geometry_SCF_conv_tol"]
        return mean_field

    def optimizer(backend):
        driver = original_optimizer(backend)

        def run(method, **kwargs):
            result = driver(method, **kwargs)
            checkpoint_directory.mkdir(parents=True, exist_ok=True)
            write_geometry(
                result[1],
                checkpoint_directory / "optimizer_returned.xyz",
                f"Optimizer returned; convergence={bool(result[0])}; "
                "physical QC still required",
            )
            return result

        return run

    try:
        dft.dft_smd = evaluate
        dft._build_smd_mf = build
        dft._geomopt_fn = optimizer
        yield
    finally:
        dft.dft_smd = original_dft
        dft._build_smd_mf = original_builder
        dft._geomopt_fn = original_optimizer


def calculate_energy(campaign, task, pilot=False):
    campaign = Path(campaign)
    plan = read_plan(campaign)
    backend = load_module(
        "production_oligomer_backend", ROOT / "scripts/pipeline/run_pi_stacking.py"
    )
    backend.BASE = campaign / "energy"
    energy_plan = dict(plan, raw_root=str(Path(plan["raw_root"]) / "energy"))
    checkpoint = (
        Path(energy_plan["raw_root"]) / task["key"] / "rest/geometry_checkpoints"
    )
    with optimizer_settings(plan["quantum_settings"], checkpoint):
        backend.calculate(task, energy_plan, pilot=pilot)
    if pilot:
        return
    audit = load_module(
        "polymer_physical_audit",
        ROOT / "scripts/validation/audit_polymer_geometries.py",
    )
    review = audit.audit(
        task, Path(energy_plan["raw_root"]) / task["key"] / "rest/opt.xyz"
    )
    if review["issues"]:
        raise ValueError(f"Physical geometry QC failed: {review['issues']}")
    directory = campaign / "energy/tasks" / task["key"]
    atomic_write_json(directory / "geometry_review.json", review)
    report_path = directory / "status.json"
    report = json.loads(report_path.read_text())
    report["output_sha256"][str(directory / "geometry_review.json")] = digest(
        directory / "geometry_review.json"
    )
    report["physical_QC"] = "passed"
    atomic_write_json(report_path, report)


def calculate_native_G(campaign, energy_record):
    """Use separate state storage so energy completion evidence stays immutable."""
    from redox.qm.sp import ensure
    import redox.qm.native_thermal as native
    from redox.core.protocol import load_record
    from redox.core.common import apply_thermal_policy, free_energy
    from redox.workflows.resources import verify_gate_owner

    verify_gate_owner()
    campaign = Path(campaign)
    plan = read_plan(campaign)
    task = energy_record["task"]
    key = task["key"]
    directory = campaign / "thermal/tasks" / key
    directory.mkdir(exist_ok=False)
    report = dict(
        task=task,
        status="computing",
        stage="native_thermal",
        started=time.time(),
        parent_energy_state=str(campaign / "energy/tasks" / key / "status.json"),
    )
    atomic_write_json(directory / "status.json", report)
    state = Path(plan["raw_root"]) / "thermal" / key / "rest"
    state.mkdir(parents=True, exist_ok=False)
    source = Path(energy_record["state_dir"])
    for name in ("opt.xyz", "result.json"):
        shutil.copyfile(source / name, state / name)
    ensure(key, "rest", state, 0, 1, backend="gpu")
    manifest_row = dict(
        id=key, state="rest", charge=0, mult_hint=1, smiles=task["smiles"]
    )
    previous_manifest = native.read_manifest
    try:
        native.read_manifest = lambda: [manifest_row]
        native.select_policy(state)
        native.ensure_native_thermal(
            key,
            "rest",
            state_dir=state,
            memory_cap_gib=plan["thermochemistry"]["memory_cap_gib"],
        )
    finally:
        native.read_manifest = previous_manifest
    raw = json.loads((state / "result.json").read_text())
    active = load_record(state, 0, 1)
    if active is None:
        raise ValueError("Native thermochemistry lacks an admitted electronic energy")
    qualified = apply_thermal_policy(dict(raw, e_smd_eV=active["e_smd_eV"]), state)
    value, qc = free_energy(qualified)
    if value is None or qc != "ok":
        raise ValueError(f"Native free-energy QC failed: {qc}")
    verify_sources(plan["producer_sha256"])
    outputs = {
        str(path): digest(path)
        for path in state.rglob("*")
        if path.is_file() and path.suffix in (".json", ".xyz", ".npy")
    }
    report.update(
        status="state_records_verified",
        G_eV=value,
        E_smd_eV=active["e_smd_eV"],
        state_dir=str(state),
        thermal_QC=qc,
        thermal_definition=qualified["thermal_definition"],
        output_sha256=outputs,
        finished=time.time(),
    )
    atomic_write_json(directory / "status.json", report)
