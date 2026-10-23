#!/usr/bin/env python
"""Check exact enumeration/route identity and every prepared polymer geometry."""
import argparse
import json
import os
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import rdDetermineBonds

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.build.attachment_identity import exact_attachment
from redox.build.grafting import precursor
from redox.build.geometry import _skeleton
from redox.build.polymer_contacts import pi_core_groups
from redox.workflows import polymer_campaign as api
from redox.workflows.resources import allocation_guard, cpu_allocation, host_capacity


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--campaign", type=Path, required=True)
    args = parser.parse_args()
    with allocation_guard():
        host_capacity(1)
        os.sched_setaffinity(0, cpu_allocation(1))
    base = args.campaign.resolve()
    plan = json.loads((base / "plan.json").read_text())
    api.verify_sources(plan["source_sha256"])
    candidate = json.loads((base / "candidate_snapshot.json").read_text())
    route = json.loads((base / "route_snapshot.json").read_text())
    fragment = candidate["frag"]
    pre, _ = precursor(fragment)
    states = candidate["states"]
    states = json.loads(states) if isinstance(states, str) else states
    binding = exact_attachment(dict(identity=dict(attachment_fragment=fragment,
        scaffold_smiles="Cc1ccc(C[*:1])cc1", product_smiles=candidate["exact_grafted_smiles"], precursor_smiles=Chem.MolToSmiles(pre)),
        payload=dict(attachment_fragment=fragment, product_smiles=candidate["exact_grafted_smiles"],
            precursor_smiles=Chem.MolToSmiles(pre), family=candidate["family"], states=states)))["route_binding"]
    if binding != route["candidate"]["route_binding"] or not route["gate"]["proceed_to_feasibility_model"]:
        raise ValueError("Route evidence does not bind this exact precursor/product/new bond")
    audit_source = ROOT / "scripts/validation/audit_polymer_geometries.py"
    audit = api.load_module("candidate_occupancy_geometry_review", audit_source)
    reviews = []
    sources = dict(plan["source_sha256"])
    for task in plan["tasks"]:
        case_path = Path(task["case"])
        xyz = Path(task["input"])
        case = json.loads(case_path.read_text())
        observed = Chem.MolFromXYZFile(str(xyz))
        rdDetermineBonds.DetermineConnectivity(observed, charge=0)
        declared = Chem.AddHs(Chem.MolFromSmiles(case["smiles"]))
        mapping = _skeleton(observed).GetSubstructMatch(_skeleton(declared))
        if len(mapping) != declared.GetNumAtoms():
            raise ValueError("Prepared XYZ does not map to the declared full atom graph")
        groups = [[mapping[i] for i in group] for group in pi_core_groups(declared, fragment)]
        review = audit.audit(dict(case, model_type="covalent_oligomer", pi_groups=groups), xyz)
        if review["issues"] or len(groups) != len(case["occupied_sites"]):
            raise ValueError("Prepared candidate geometry failed physical/topology QC")
        if (case["backbone_carbon_count"], case["backbone_CC_bonds"], case["chain_ends"]) != (10, 9, 2):
            raise ValueError("Unexpected polymer backbone topology")
        reviews.append(review)
        sources.update({str(p): api.digest(p) for p in (case_path, xyz)})
    sources.update({str(p): api.digest(p) for p in (Path(__file__).resolve(), base / "plan.json", audit_source)})
    atomic_write_json(base / "input_validation.json", dict(status="exact_candidate_occupancy_inputs_checked",
        candidate_id=candidate["id"], enumerated_aliases=plan["enumerated_aliases"], route_binding=binding,
        precursor_access_reviewed=True, prospective_route_status=plan["prospective_route_status"],
        geometry_reviews=reviews, source_sha256=sources, full_redox_admission=False, synthesis_admission=False,
        minimum_Hessian_proven=False, energies_computed=False))
    print(json.dumps(dict(candidate=candidate["id"], exact_pool_aliases=plan["enumerated_aliases"], checked_geometries=len(reviews), route_status=plan["prospective_route_status"])))


if __name__ == "__main__":
    main()
