#!/usr/bin/env python
"""Freeze checked candidate MMFF starts for an explicitly limited DFT baseline.

These are independent warm starts, not intermediate MD frames or a converged
conformer ensemble. The broader source searches continue. Case directories share
the existing resource root so every launcher accounts for the same RAM promises.
"""
import argparse
import json
import os
from pathlib import Path
import shutil
import time

from ase.io import read
from rdkit import Chem
from rdkit.Chem import rdDetermineBonds
from pyscf import gto

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.core.protocol import ACTIVE_SP
from redox.build.geometry import _skeleton
from redox.build.polymer_contacts import pi_core_groups
from redox.workflows import polymer_campaign as api
from redox.workflows.resources import allocation_guard, cpu_allocation, host_capacity


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--campaign", type=Path, required=True)
    parser.add_argument("--resource-root", type=Path, required=True)
    args = parser.parse_args()
    with allocation_guard():
        host_capacity(1)
        os.sched_setaffinity(0, cpu_allocation(1))
    base = args.campaign.resolve()
    resource_root = args.resource_root.resolve()
    if (base / "baseline_plan.json").exists():
        raise FileExistsError("Preserve the existing baseline declaration")
    inputs = json.loads((base / "plan.json").read_text())
    review = json.loads((base / "input_validation.json").read_text())
    if review["status"] != "exact_candidate_occupancy_inputs_checked":
        raise ValueError("Candidate inputs need exact identity/route/geometry review")
    api.verify_sources(review["source_sha256"])
    candidate = json.loads((base / "candidate_snapshot.json").read_text())
    parent = ROOT / "results/polymer_models/campaigns/2026-10-21_spacing_cds_recovery"
    sources = list((ROOT / "src/redox").rglob("*.py")) + [Path(__file__).resolve()]
    sources += [ROOT / "scripts/pipeline" / name for name in (
        "refine_polymer_conformers.py", "run_polymer_spacing.py", "run_pi_stacking.py", "polymer_conformer_thermals.py",
        "resume_polymer_forces_adaptive.py", "resume_polymer_forces.py", "continue_polymer_models.py",
        "resume_bounded_polymer_thermal.py", "resume_polymer_thermal.py")]
    sources += [ROOT / "scripts/validation" / name for name in ("audit_polymer_geometries.py", "diagnose_saved_hessian.py")]
    cases = []
    audit = api.load_module("candidate_baseline_geometry", ROOT / "scripts/validation/audit_polymer_geometries.py")
    for prepared in inputs["tasks"]:
        case = json.loads(Path(prepared["case"]).read_text())
        key = "baseline__" + prepared["key"]
        campaign = resource_root / key
        campaign.mkdir(exist_ok=False)
        task_key = key + "__pose0"
        directory = campaign / "live_inputs" / task_key
        directory.mkdir(parents=True)
        xyz = directory / "input.xyz"
        shutil.copyfile(prepared["input"], xyz)
        atoms = read(xyz)
        declared = Chem.AddHs(Chem.MolFromSmiles(case["smiles"]))
        observed = Chem.MolFromXYZFile(str(xyz))
        rdDetermineBonds.DetermineConnectivity(observed, charge=0)
        mapping = _skeleton(observed).GetSubstructMatch(_skeleton(declared))
        if len(mapping) != declared.GetNumAtoms():
            raise ValueError("Candidate baseline input atom mapping failed")
        groups = [[mapping[i] for i in group] for group in pi_core_groups(declared, candidate["frag"])]
        basis = gto.M(atom=list(zip(atoms.get_chemical_symbols(), atoms.positions)), basis=ACTIVE_SP["basis"], charge=0, spin=0, unit="Angstrom", verbose=0).nao_nr()
        task = dict(case, key=task_key, model_key=key, candidate_id=candidate["id"], model_type="covalent_oligomer",
            repeats=case["repeat_units"], stereoisomer="all_same", target="validated_MMFF_start_baseline",
            fragment=candidate["frag"], pi_groups=groups, active_basis_functions=basis,
            search_rank=1, source_conformer_index=None, source_GFN2_energy_Ha=None,
            starting_geometry_source_role="validated_MMFF_baseline; broad_CREST_concurrent")
        geometry = audit.audit(task, xyz)
        if geometry["issues"]:
            raise ValueError("Prepared baseline geometry failed physical input QC")
        atomic_write_json(directory / "geometry_review.json", geometry)
        hashes = {str(path): api.digest(path) for path in (xyz, Path(prepared["input"]), Path(prepared["case"]), base / "input_validation.json", directory / "geometry_review.json")}
        atomic_write_json(directory / "task.json", dict(task, input_xyz=str(xyz), input_sha256=hashes))
        hashes[str(directory / "task.json")] = api.digest(directory / "task.json")
        for stage in ("energy", "thermal"):
            for name in ("tasks", "done", "claims", "failures", "workers", "pilot", "launches"):
                (campaign / stage / name).mkdir(parents=True)
        shutil.copytree(parent / "controls", campaign / "controls")
        for path in (campaign / "controls").glob("*/status.json"):
            hashes[str(path)] = api.digest(path)
        plan = dict(schema=1, tasks=[task], graph_cases=[dict(case, key=key)], cycles=[],
            raw_root=str(ROOT / "calcs/polymer_models" / base.name / "baseline" / key),
            source_search=prepared["raw_out"], source_role="checked_MMFF_warm_start; not a ranked CREST minimum",
            ordered_DFT_batches=[10, 20], quantum_settings=dict(max_opt_steps=400, geometry_SCF_conv_tol=1e-11),
            thermochemistry=dict(temperature_K=298.15, memory_cap_gib=24), solvent="acetonitrile",
            minimum_large_GPU_MiB=80 * 1024, maximum_memory_pilot_MiB=80 * 1024, pilot_key=task_key,
            producer_sha256={str(path): api.digest(path) for path in sources}, input_sha256=hashes,
            created=time.time(), synthesis_admission=False, conformer_convergence_claimed=False)
        atomic_write_json(campaign / "plan.json", plan)
        cases.append(dict(key=key, campaign=str(campaign), occupied_sites=case["occupied_sites"], active_basis_functions=basis, atoms=case["atoms"]))
    atomic_write_json(base / "cases.json", dict(fragment=candidate["frag"], motif_id=candidate["id"]))
    terms = api.reaction_terms(base, candidate["id"], "dbu")
    atomic_write_json(base / "baseline_plan.json", dict(cases=cases, largest_case=max(cases, key=lambda row: row["active_basis_functions"])["key"],
        resource_root=str(resource_root), shared_empty_reference=inputs["common_empty_reference"],
        reaction_terms_status=terms["status"], scope="One-pose occupied-case feasibility baseline; wider conformer ensemble remains pending",
        source_sha256={str(path): api.digest(path) for path in (Path(__file__).resolve(), base / "plan.json", base / "input_validation.json", base / "reaction_terms.json")}))
    print(json.dumps(dict(cases=len(cases),largest=max(cases,key=lambda row:row["active_basis_functions"]))))


if __name__ == "__main__":
    main()
