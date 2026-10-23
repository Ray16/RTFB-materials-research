#!/usr/bin/env python
"""Recompute segment deformation and inspect saved optimization/geometry evidence.

Read-only molecular audit; no new QM or relaxation. A small deformation energy
cannot establish a conformational ensemble or isolate a pendant-core interaction.
"""
import argparse
import json
import os
from pathlib import Path
import re
import time

import numpy as np
from ase.io import read
from rdkit import Chem
from rdkit.Chem import rdMolAlign

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.properties.conformer_ensemble import EV_KJ_MOL
from redox.workflows import polymer_campaign as api
from redox.workflows.checkpoints import verify_state_checkpoint
from redox.workflows.resources import allocation_guard, cpu_allocation, host_capacity


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    with allocation_guard():
        host_capacity(1)
        os.sched_setaffinity(0, cpu_allocation(1))
    args.out.mkdir(parents=True, exist_ok=False)
    geometry = api.load_module("segment_geometry_mapping", ROOT / "scripts/pipeline/polymer_conformer_thermals.py")
    fragments = ROOT / "results/polymer_models/campaigns/2026-10-20_frozen_fragments"
    original = ROOT / "results/polymer_models/campaigns/2026-10-19_pi_stacking"
    decomposition = ROOT / "results/polymer_models/interaction_decomposition.json"
    checked = ROOT / "results/polymer_models/interaction_checked.json"
    source_paths = {Path(__file__).resolve(), decomposition, checked, fragments / "plan.json", ROOT / "scripts/pipeline/polymer_conformer_thermals.py"}
    models = {row["key"]: row for row in json.loads(checked.read_text())["models"]}
    tasks = {row["key"]: row for row in api.read_plan(fragments)["tasks"]}

    def parent(row):
        campaign = Path(row["recovery_campaign"]) if row["recovery_campaign"] else original
        if not campaign.is_absolute():
            campaign = ROOT / campaign
        receipt = campaign / "done" / (row["key"] + ".json")
        report = verify_state_checkpoint(receipt)
        source_paths.add(receipt)
        if abs(report["E_smd_eV"] - row["E_smd_eV"]) > 1e-8:
            raise ValueError("Published segment energy differs from its receipt")
        xyz = Path(row["final_xyz"])
        raw = json.loads((xyz.parent / "result.json").read_text())
        if not all(raw.get(name) is True for name in ("geometry_optimization_converged", "optimization_SCF_converged", "opt_gradient_grid_response")):
            raise ValueError("Relaxed reference lacks optimizer/SCF/full-grid convergence")
        source_paths.update((xyz, xyz.parent / "result.json", Path(report["task"]["input_xyz"])))
        log = campaign / "tasks" / row["key"] / "calculation.log"
        steps = None
        if log.exists():
            source_paths.add(log)
            text = re.sub(r"\x1b\[[0-9;]*m", "", log.read_text(errors="replace"))
            matches = [int(value) for value in re.findall(r"^Step\s+(\d+)\s*:", text, flags=re.MULTILINE)]
            steps = max(matches) if matches else None
        return report, dict(optimizer_steps_in_saved_log=steps, gradient_max_Ha_Bohr=report["gradient_max_Ha_Bohr"], Hessian_minimum_tested=False)

    pairs = []
    for pair in json.loads(decomposition.read_text())["pairs"]:
        model = models[pair["pair_key"]]
        pair_parent, pair_optimization = parent(model)
        monomer = next(row for row in models.values() if row["candidate_id"] == pair["candidate_id"] and row["target"] == "monomer")
        mono_parent, mono_optimization = parent(monomer)
        mono_xyz = Path(monomer["final_xyz"])
        reference = geometry.alignment_molecule(monomer["smiles"], mono_xyz)
        pair_xyz = read(model["final_xyz"])
        frozen_rows = []
        for key in pair["fragment_keys"]:
            receipt = fragments / "done" / (key + ".json")
            report = verify_state_checkpoint(receipt)
            task = tasks[key]
            xyz = Path(task["input_xyz"])
            frozen = read(xyz)
            assert frozen.get_chemical_symbols() == pair_xyz[task["source_atom_indices"]].get_chemical_symbols()
            if not np.allclose(frozen.positions, pair_xyz.positions[task["source_atom_indices"]], atol=1e-8, rtol=0):
                raise ValueError("Frozen fragment does not match the extracted optimized pair")
            molecule = geometry.alignment_molecule(task["smiles"], xyz)
            rmsd = rdMolAlign.GetBestRMS(Chem.Mol(molecule), reference, maxMatches=10000)
            frozen_rows.append(dict(key=key, atoms=task["atoms"], smiles=task["smiles"],
                E_eV=report["E_smd_eV"], deformation_kJ_mol=(report["E_smd_eV"] - mono_parent["E_smd_eV"]) * EV_KJ_MOL,
                frozen_to_relaxed_monomer_heavy_RMSD_A=rmsd))
            source_paths.update((receipt, xyz))
        pair_initial = geometry.alignment_molecule(model["smiles"], pair_parent["task"]["input_xyz"])
        pair_final = geometry.alignment_molecule(model["smiles"], model["final_xyz"])
        relaxation_rmsd = rdMolAlign.GetBestRMS(pair_initial, pair_final, maxMatches=10000)
        total = sum(row["deformation_kJ_mol"] for row in frozen_rows)
        stored_constant_replay = total * 96.485 / EV_KJ_MOL
        if abs(stored_constant_replay - pair["segment_deformation_deltaE_kJ_mol"]) > 1e-6:
            raise ValueError("Deformation energy does not reconstruct")
        pairs.append(dict(pair_key=pair["pair_key"], candidate_id=pair["candidate_id"],
            pair_atoms=model["atoms"], pair_optimization=pair_optimization, monomer_optimization=mono_optimization,
            pair_input_to_optimized_heavy_RMSD_A=relaxation_rmsd, fragments=frozen_rows,
            deformation_kJ_mol_CODATA=total, historical_deformation_kJ_mol=pair["segment_deformation_deltaE_kJ_mol"],
            frozen_interaction_kJ_mol=pair["frozen_interaction_deltaE_kJ_mol"],
            negative_fragment_terms=[row["key"] for row in frozen_rows if row["deformation_kJ_mol"] < 0],
            cancellation_kJ_mol=sum(abs(row["deformation_kJ_mol"]) for row in frozen_rows) - abs(total),
            deformation_is_positive_definite_strain=False,
            interpretation="Signed difference from one sampled independently optimized monomer. Small sums may include cancellation; geometry and single-point levels differ. Neither global monomer minima nor nonnegative strain penalties are established."))
    source = {str(path): api.digest(path) for path in source_paths}
    atomic_write_json(args.out / "summary.json", dict(status="segment_deformation_and_saved_relaxation_audited", checked=time.time(), pairs=pairs,
        source_sha256=source, new_quantum_calculations=0, source_data_modified=False,
        original_constraints_reconstructed=False,
        limitations="Stationary electronic structures only; no new Hessian, global minimum, constraint-history reconstruction, BSSE correction or pure-core binding result"))
    (args.out / "producer_source.py").write_bytes(Path(__file__).read_bytes())
    print(json.dumps(dict(pairs=len(pairs), deformation_kJ_mol=[row["deformation_kJ_mol_CODATA"] for row in pairs],
        pair_relaxation_RMSD_A=[row["pair_input_to_optimized_heavy_RMSD_A"] for row in pairs])))


if __name__ == "__main__":
    main()
