"""Export/replay SI controls while preserving their distinct scientific roles.

Only preparation imports the scientific environment. Cached arithmetic uses the
standard library; alternative-basis records never acquire a primary certificate.
"""
import json
from pathlib import Path


def validate_alternative_surface(surface):
    if (surface.get("primary_native_certificate_claimed") is not False
            or surface["identity"]["protocol"].get("basis") != "def2-svpd"
            or surface.get("status") != "derivative_checked_solution_thermal"
            or not surface.get("stationary_minimum") or not surface.get("grid_response")
            or not 0 <= surface["gradient_max"] <= .00045
            or not 0 <= surface["maximum_asymmetry"] <= 1e-6):
        raise ValueError("Alternative basis scope or original physical gates differ")
    if len(surface["derivative_checks"]) < 4 or any(r["relative_error"] > .01 for r in surface["derivative_checks"]):
        raise ValueError("Alternative basis derivative evidence is incomplete")


def add_auxiliary_states(root, load, state, sources, digest):
    from redox.workflows import polymer_campaign as api
    from redox.core.common import read_manifest
    manifest = {(r["id"], r["state"]): r for r in read_manifest()}
    result = dict(proof_sources=[], basis=[], restart=[])
    basis = load("neutral_basis.json")
    api.verify_sources(basis["source_sha256"])
    for row in basis["rows"]:
        surface = Path(row["alternative_surface"])
        status = json.loads((surface.parent / "status.json").read_text())
        if status["status"] != "alternative_basis_control_checked" or status["primary_native_admission"] is not False:
            raise ValueError("Unqualified or relabeled basis control")
        validate_alternative_surface(json.loads(surface.read_text()))
        state("basis:" + row["id"], surface.parent, manifest[row["id"], "neu"]["smiles"],
            status["G_eV"], alternative_surface=surface)
        result["basis"].append(dict(id=row["id"], alias="basis:" + row["id"], primary_alias="benchmark:" + row["id"] + ":neu"))
        result["proof_sources"].append(str(surface.parent / "status.json"))
    campaign = root / "results/validation/quantum_restart_2026-10-22"
    plan_path = campaign / "plan.json"
    plan = api.read_plan(campaign)
    if digest(plan_path) != load("quantum_restart.json")["plan_sha256"]:
        raise ValueError("Restart plan differs from the frozen control")
    for task in plan["tasks"]:
        folder = Path(plan["raw_root"]) / task["key"]
        status_path = folder / "status.json"
        status = json.loads(status_path.read_text())
        result["proof_sources"].append(str(status_path))
        if task["phase"] == "interrupted":
            proof_path = campaign / "death_checks" / (task["key"] + ".json")
            proof = json.loads(proof_path.read_text())
            if not proof["owner_dead"] or proof["expected_interrupt_returncode"] != 87:
                raise ValueError("Saved interruption ownership proof differs")
            checkpoint = folder / "geometry_checkpoints/last_evaluated.xyz"
            if digest(checkpoint) != proof["checkpoint_sha256"]:
                raise ValueError("Saved unqualified checkpoint differs")
            result["proof_sources"] += [str(proof_path), str(checkpoint)]
            continue
        if status["status"] != "quantum_control_verified":
            raise ValueError("Restart control has not completed")
        api.verify_sources(status["output_sha256"])
        state("restart:" + task["key"], status["state_dir"], task["smiles"], status["G_eV"])
        result["restart"].append(dict(case=task["case"], phase=task["phase"], alias="restart:" + task["key"]))
    for path in result["proof_sources"]:
        sources[path] = digest(path)
    return result


def electronic_inventory(root, load, artifact, digest):
    from redox.core.protocol import load_record, record_path
    from redox.workflows.checkpoints import verify_state_checkpoint
    current_path = root / "results/polymer_models/interaction_checked.json"
    current = json.loads(current_path.read_text())
    accepted = [r for r in current["models"] if r["stage"] == "stationary_DFT_electronic_structure" and r["E_smd_eV"] is not None]
    if len(accepted) != load("polymer_validation.json")["original_checked"]:
        raise ValueError("Original electronic model coverage differs")
    original = root / "results/polymer_models/campaigns/2026-10-19_pi_stacking"
    fragments = root / "results/polymer_models/campaigns/2026-10-20_frozen_fragments"
    rows = []
    def record(key, receipt, geometry, expected, role):
        report = verify_state_checkpoint(receipt)
        sd = Path(report["state_dir"])
        raw = json.loads((sd / "result.json").read_text())
        active = load_record(sd, raw["charge"], raw["mult"])
        if active is None or abs(active["e_smd_eV"] - expected) > 1e-8 or abs(report["E_smd_eV"] - expected) > 1e-8:
            raise ValueError("Electronic-only source energy differs")
        if digest(geometry) != digest(sd / "opt.xyz"):
            raise ValueError("Electronic-only source geometry differs")
        receipt_data = json.loads(receipt.read_text())
        paths = dict(geometry=geometry, result=sd / "result.json", energy=record_path(sd, raw["charge"], raw["mult"]),
            receipt=receipt, status=Path(receipt_data["status_path"]), initial_geometry=Path(report["task"]["input_xyz"]))
        rows.append(dict(id=key, role=role, E_eV=expected, charge=raw["charge"], mult=raw["mult"],
            files={name:artifact(path) for name,path in paths.items()}, free_energy_available=False,
            Hessian_minimum_established=False))
    for row in accepted:
        campaign = Path(row["recovery_campaign"]) if row["recovery_campaign"] else original
        if not campaign.is_absolute():
            campaign = root / campaign
        record(row["key"], campaign / "done" / (row["key"] + ".json"), Path(row["final_xyz"]),
            row["E_smd_eV"], "stationary_electronic_model")
    for pair in load("segment_relaxation.json")["pairs"]:
        for fragment in pair["fragments"]:
            key = fragment["key"]
            receipt = fragments / "done" / (key + ".json")
            report = verify_state_checkpoint(receipt)
            record(key, receipt, Path(report["task"]["input_xyz"]), fragment["E_eV"], "frozen_fragment_single_point")
    return dict(records=rows, checked_models=artifact(current_path), excluded_models=[
        dict(id=r["key"], status=r["stage"], E_eV=None, free_energy_available=False)
        for r in current["models"] if r not in accepted])


def replay_auxiliary(root, data, g, load, local):
    """Arithmetic on archived inputs; no original path or quantum import."""
    aux = data["auxiliary"]
    rows = {r["id"]:r for r in data["states"]}
    aliases = data["aliases"]
    read = lambda relative: json.loads(local(root / "data", relative).read_text())
    energy = lambda alias: read(rows[aliases[alias]]["files"]["energy"])
    thermal = lambda alias: read(rows[aliases[alias]]["files"]["partition"])
    basis = {r["id"]:r for r in load("neutral_basis.json")["rows"]}
    for row in aux["basis"]:
        actual, primary = row["alias"], row["primary_alias"]
        surface = read(rows[aliases[actual]]["files"]["native"])
        validate_alternative_surface(surface)
        shifts = dict(delta_thermal_eV=thermal(actual)["g_thermal_total_eV"]-thermal(primary)["g_thermal_total_eV"],
            delta_geometry_SP_eV=energy(actual)["e_smd_eV"]-energy(primary)["e_smd_eV"],
            delta_first_reduction_V_holding_anion_and_Fc_fixed=g(actual)-g(primary))
        if any(abs(value-basis[row["id"]][key]) > 1e-8 for key,value in shifts.items()):
            raise ValueError("Portable basis-control shifts differ")
    for control in load("quantum_restart.json")["controls"]:
        prefix = "restart:" + control["case"]
        a, b = prefix + "__baseline", prefix + "__resumed"
        ea, eb = energy(a), energy(b)
        difference = max(abs(ea[k]-eb[k]) for k in ("e_smd_eV", "e_gas_eV"))
        if (abs(difference-control["differences"]["energy_eV"]) > 1e-9
                or abs(abs(g(a)-g(b))-control["differences"]["free_energy_eV"]) > 1e-9):
            raise ValueError("Portable quantum-restart energy comparison differs")
    electronic = {r["id"]:r for r in aux["electronic"]["records"]}
    geometry = {}
    for row in electronic.values():
        active = read(row["files"]["energy"])
        raw = read(row["files"]["result"])
        blob = local(root / "data", row["files"]["geometry"]).read_bytes()
        import hashlib
        if (not active["converged_smd"] or hashlib.sha1(blob).hexdigest()[:16] != active["geom_sha1"]
                or abs(active["e_smd_eV"]-row["E_eV"]) > 1e-8 or row["free_energy_available"]
                or any(value.get("charge") != row["charge"] or value.get("mult") != row["mult"] for value in (raw, active))):
            raise ValueError("Electronic-only geometry/energy identity or scope differs")
        geometry[row["id"]] = blob.decode().splitlines()[2:]
    models = read(aux["electronic"]["checked_models"])["models"]
    for pair in load("segment_decomposition.json")["pairs"]:
        complex_energy = electronic[pair["pair_key"]]["E_eV"]
        model = next(r for r in models if r["candidate_id"] == pair["candidate_id"] and r["target"] == "monomer")
        mono = electronic[model["key"]]["E_eV"]
        fragments = sum(electronic[key]["E_eV"] for key in pair["fragment_keys"])
        checks = ((complex_energy-fragments)*96.485, (fragments-2*mono)*96.485, (complex_energy-2*mono)*96.485)
        names = ("frozen_interaction_deltaE_kJ_mol", "segment_deformation_deltaE_kJ_mol", "relaxed_association_deltaE_kJ_mol")
        if any(abs(value-pair[name]) > 1e-6 for name,value in zip(names,checks)):
            raise ValueError("Portable electronic decomposition differs")
        for key in pair["fragment_keys"]:
            status = read(electronic[key]["files"]["status"])
            indexes = status["task"]["source_atom_indices"]
            a = [line.split() for line in geometry[key]]
            b = [geometry[pair["pair_key"]][i].split() for i in indexes]
            if len(a) != len(b) or any(x[0]!=y[0] or max(abs(float(v)-float(w)) for v,w in zip(x[1:],y[1:])) > 1e-8 for x,y in zip(a,b)):
                raise ValueError("Frozen fragment does not match archived parent coordinates")
    return dict(basis_controls=len(aux["basis"]), restart_controls=2,
        stationary_electronic_models=sum(r["role"]=="stationary_electronic_model" for r in electronic.values()),
        frozen_fragment_single_points=sum(r["role"]=="frozen_fragment_single_point" for r in electronic.values()),
        electronic_decomposition_pairs=len(load("segment_decomposition.json")["pairs"]))
